//! `HtsgetReader` against the routes the dump served, on the same two ports: TLS on 18443 and plain
//! on 18080. One test, since the working directory and the trusted roots are the process's.
//!
//! The golden is absent until a CI artefact freezes it, and a golden that is absent passes.

use gatk_corpus as corpus;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

struct Route {
    statuses: Vec<u16>,
    location: String,
    body: Vec<u8>,
    seen: AtomicUsize,
}

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// One request off the stream, and its answer. The path, and the two headers `/echo` reports.
fn answer<S: std::io::Read + Write>(stream: S, routes: &HashMap<String, Route>) {
    let mut reader = BufReader::new(stream);
    let mut line = String::new();
    if reader.read_line(&mut line).is_err() {
        return;
    }
    let target = line.split(' ').nth(1).unwrap_or("/").to_string();
    let path = target.split('?').next().unwrap_or("/").to_string();
    let (mut test, mut number) = ("null".to_string(), "null".to_string());
    loop {
        line.clear();
        if reader.read_line(&mut line).unwrap_or(0) == 0 || line == "\r\n" {
            break;
        }
        let Some((name, value)) = line.split_once(':') else {
            continue;
        };
        match name.to_ascii_lowercase().as_str() {
            "x-test" => test = value.trim().to_string(),
            "x-number" => number = value.trim().to_string(),
            _ => {}
        }
    }
    let (status, location, body) = match routes.get(&path) {
        Some(route) => {
            let seen = route.seen.fetch_add(1, Ordering::SeqCst);
            let status = route.statuses[seen.min(route.statuses.len() - 1)];
            let body = if path == "/echo" {
                format!("{test}/{number}\n").into_bytes()
            } else {
                route.body.clone()
            };
            (status, route.location.clone(), body)
        }
        None => (404, "-".to_string(), b"no route".to_vec()),
    };
    let mut head = format!(
        "HTTP/1.1 {status} Answer\r\nContent-Length: {}\r\nConnection: close\r\n",
        body.len()
    );
    if location != "-" {
        head.push_str(&format!("Location: {location}\r\n"));
    }
    head.push_str("\r\n");
    let stream = reader.get_mut();
    let _ = stream.write_all(head.as_bytes());
    let _ = stream.write_all(&body);
    let _ = stream.flush();
}

#[test]
fn every_run_replays_as_the_reference_ran_it() {
    let golden = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../gatk-tools/tests/data/htsget_reader.txt.gz");
    if !golden.exists() {
        return;
    }
    let text = corpus::read_golden(&golden);

    let root = std::env::temp_dir().join("gatk-cli-htsget");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(root.join("htsget-dump/out")).expect("a scratch directory");
    std::env::set_current_dir(&root).expect("the scratch directory");

    let mut routes = HashMap::new();
    let mut labels = Vec::new();
    for line in text.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[..] {
            ["route", path, statuses, location, ref body @ ..] => {
                routes.insert(
                    path.to_string(),
                    Route {
                        statuses: statuses.split(',').map(|s| s.parse().unwrap()).collect(),
                        location: location.to_string(),
                        body: corpus::decode_base64(body.first().copied().unwrap_or("")),
                        seen: AtomicUsize::new(0),
                    },
                );
            }
            ["argv", label, ..] => labels.push(label.to_string()),
            _ => {}
        }
    }
    let routes = Arc::new(routes);

    // The dump makes its certificate with keytool; this one is made here, and trusted the same way.
    let certified = rcgen::generate_simple_self_signed(vec!["localhost".to_string()]).unwrap();
    std::fs::write("roots.pem", certified.cert.pem()).unwrap();
    std::env::set_var("SSL_CERT_FILE", root.join("roots.pem"));
    let config = Arc::new(
        rustls::ServerConfig::builder()
            .with_no_client_auth()
            .with_single_cert(
                vec![certified.cert.der().clone()],
                rustls::pki_types::PrivateKeyDer::Pkcs8(certified.key_pair.serialize_der().into()),
            )
            .unwrap(),
    );
    let secure = std::net::TcpListener::bind(("localhost", 18443)).expect("the TLS port");
    let plain = std::net::TcpListener::bind(("localhost", 18080)).expect("the plain port");
    {
        let routes = Arc::clone(&routes);
        std::thread::spawn(move || {
            for stream in plain.incoming().flatten() {
                let routes = Arc::clone(&routes);
                std::thread::spawn(move || answer(stream, &routes));
            }
        });
    }
    {
        let routes = Arc::clone(&routes);
        std::thread::spawn(move || {
            for stream in secure.incoming().flatten() {
                let routes = Arc::clone(&routes);
                let config = Arc::clone(&config);
                std::thread::spawn(move || {
                    let connection = rustls::ServerConnection::new(config).unwrap();
                    answer(rustls::StreamOwned::new(connection, stream), &routes);
                });
            }
        });
    }

    let row = |kind: &str, label: &str| -> Option<String> {
        let prefix = format!("{kind}\t{label}\t");
        text.lines()
            .find(|line| line.starts_with(&prefix))
            .map(|line| line[prefix.len()..].to_string())
    };
    let failures = Mutex::new(Vec::new());
    assert!(!labels.is_empty(), "the golden holds no runs");
    for label in &labels {
        let argv: Vec<String> = unescape(&row("argv", label).expect("a command line"))
            .split('\t')
            .map(str::to_string)
            .collect();
        let ran = gatk_cli::run(&argv);
        let mut failures = failures.lock().unwrap();
        match row("error", label) {
            Some(error) => {
                let error = unescape(&error);
                let message = match error.split_once(':') {
                    Some((_, "null")) | None => "",
                    Some((_, message)) => message,
                };
                if ran.status == 0 || !ran.stderr.contains(message) {
                    failures.push(format!(
                        "{label}: status {}, reference refused with {error:?}; stderr {:?}",
                        ran.status, ran.stderr
                    ));
                }
            }
            None if ran.status != 0 => {
                failures.push(format!("{label}: status {}: {}", ran.status, ran.stderr));
            }
            None => {}
        }
        let written = std::fs::read(format!("htsget-dump/out/{label}")).ok();
        let reference = row("file", label).map(|encoded| corpus::decode_base64(&encoded));
        if written != reference {
            failures.push(format!(
                "{label}: wrote {:?}, reference {:?}",
                written.map(|bytes| String::from_utf8_lossy(&bytes).into_owned()),
                reference.map(|bytes| String::from_utf8_lossy(&bytes).into_owned())
            ));
        }
    }
    let failures = failures.into_inner().unwrap();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
