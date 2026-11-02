//! The `.scorer.pkl` of `TrainVariantAnnotationsModel`'s default backend, both ways.
//!
//! **Writing**: `dill.dump` of the script's lambda is `pickle._Pickler` at protocol 4 over dill
//! 0.3.7's reducers, so the bytes are a function of the fit except for one thing. The lambda's
//! code object names the script's file, and the reference runs the script from a temporary file
//! with a random `long` in its name, so two runs of the reference differ there and nothing can
//! match it. The port writes the same pickle with that name zeroed, at the twenty digits the name
//! usually has, and the suite compares with the digits masked. The framing is `_Framer`'s: a
//! frame closes at the first `save` past 64 KiB, which is why the writer checks before each one.
//!
//! **Reading**: the opcodes are interpreted into plain values and nothing in the file runs. A
//! scorer written by something other than this script at scikit-learn 1.3.2 is refused, and
//! keeps its interpreter's refusal in the tool.
use crate::isolation_forest::{Forest, Node, Scorer, ScoringTree};
use crate::scorer_pickle_template as template;
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

/// The file name the code objects carry, where the reference's has its random digits.
pub const SCRIPT_NAME: &str = "/tmp/isolation-forest.00000000000000000000.py";

const FRAME_TARGET: usize = 64 * 1024;

const MARK: u8 = b'(';
const STOP: u8 = b'.';
const POP: u8 = b'0';
const BINFLOAT: u8 = b'G';
const BININT: u8 = b'J';
const BININT1: u8 = b'K';
const BININT2: u8 = b'M';
const NONE: u8 = b'N';
const REDUCE: u8 = b'R';
const BUILD: u8 = b'b';
const BINGET: u8 = b'h';
const LONG_BINGET: u8 = b'j';
const APPEND: u8 = b'a';
const APPENDS: u8 = b'e';
const SETITEMS: u8 = b'u';
const SETITEM: u8 = b's';
const TUPLE: u8 = b't';
const EMPTY_DICT: u8 = b'}';
const EMPTY_LIST: u8 = b']';
const EMPTY_TUPLE: u8 = b')';
const PROTO: u8 = 0x80;
const NEWOBJ: u8 = 0x81;
const TUPLE1: u8 = 0x85;
const TUPLE2: u8 = 0x86;
const TUPLE3: u8 = 0x87;
const NEWTRUE: u8 = 0x88;
const NEWFALSE: u8 = 0x89;
const LONG1: u8 = 0x8a;
const BINUNICODE: u8 = b'X';
const BINBYTES: u8 = b'B';
const SHORT_BINBYTES: u8 = b'C';
const SHORT_BINUNICODE: u8 = 0x8c;
const STACK_GLOBAL: u8 = 0x93;
const MEMOIZE: u8 = 0x94;
const FRAME: u8 = 0x95;
const GLOBAL: u8 = b'c';

// What the fixed runs leave in the memo, by the index dill gave them.
const CREATE_ARRAY: u8 = 100;
const NDARRAY: u8 = 103;
const B: u8 = 105;
const SIMPLE_B: u8 = 107;
const DTYPE: u8 = 91;
const LITTLE: u8 = 95;
const F8: u8 = 94;
const I8: u8 = 175;
const NODE_DTYPE: u8 = 191;
const VERSION_KEY: u8 = 114;
const VERSION: u8 = 115;
const TREE_CLASS: u8 = 169;
const REGRESSOR: u8 = 125;
const VALUES_KEY: u8 = 224;
const COUNT_KEY: u8 = 184;
const NODES_KEY: u8 = 185;
const DEPTH_KEY: u8 = 132;

/// `pickle._Framer` and the opcodes it is asked to write.
struct Writer {
    out: Vec<u8>,
    frame: Vec<u8>,
}

impl Writer {
    fn commit(&mut self, force: bool) {
        if self.frame.len() >= FRAME_TARGET || (force && !self.frame.is_empty()) {
            if self.frame.len() >= 4 {
                self.out.push(FRAME);
                self.out
                    .extend_from_slice(&(self.frame.len() as u64).to_le_bytes());
            }
            self.out.append(&mut self.frame);
        }
    }

    /// `save`'s first line.
    fn check(&mut self) {
        self.commit(false);
    }

    fn raw(&mut self, bytes: &[u8]) {
        self.frame.extend_from_slice(bytes);
    }

    fn op(&mut self, op: u8) {
        self.frame.push(op);
    }

    fn memoize(&mut self) {
        self.op(MEMOIZE);
    }

    fn get(&mut self, index: u8) {
        self.check();
        self.op(BINGET);
        self.op(index);
    }

    fn none(&mut self) {
        self.check();
        self.op(NONE);
    }

    fn falsity(&mut self) {
        self.check();
        self.op(NEWFALSE);
    }

    fn truth(&mut self) {
        self.check();
        self.op(NEWTRUE);
    }

    /// `save_long`.
    fn int(&mut self, value: i64) {
        self.check();
        if (0..=0xff).contains(&value) {
            self.op(BININT1);
            self.op(value as u8);
        } else if (0..=0xffff).contains(&value) {
            self.op(BININT2);
            self.raw(&(value as u16).to_le_bytes());
        } else if (-0x8000_0000..=0x7fff_ffff).contains(&value) {
            self.op(BININT);
            self.raw(&(value as i32).to_le_bytes());
        } else {
            // `encode_long`: as few little-endian two's-complement bytes as hold it.
            let bytes = value.to_le_bytes();
            let mut used = 8;
            while used > 1 {
                let top = bytes[used - 1];
                let next = bytes[used - 2];
                if (top == 0 && next & 0x80 == 0) || (top == 0xff && next & 0x80 != 0) {
                    used -= 1;
                } else {
                    break;
                }
            }
            self.op(LONG1);
            self.op(used as u8);
            self.raw(&bytes[..used]);
        }
    }

    fn float(&mut self, value: f64) {
        self.check();
        self.op(BINFLOAT);
        self.raw(&value.to_be_bytes());
    }

    /// `save_str`.
    fn string(&mut self, value: &str) {
        self.check();
        if value.len() <= 0xff {
            self.op(SHORT_BINUNICODE);
            self.op(value.len() as u8);
        } else {
            self.op(BINUNICODE);
            self.raw(&(value.len() as u32).to_le_bytes());
        }
        self.raw(value.as_bytes());
        self.memoize();
    }

    /// `save_bytes`. An array is short of the 64 KiB that would write it outside any frame: a
    /// tree here has at most 511 nodes, and the rest are smaller.
    fn bytes(&mut self, value: &[u8]) {
        self.check();
        if value.len() <= 0xff {
            self.op(SHORT_BINBYTES);
            self.op(value.len() as u8);
        } else {
            self.op(BINBYTES);
            self.raw(&(value.len() as u32).to_le_bytes());
        }
        self.raw(value);
        self.memoize();
    }

    /// `save_ndarray` in dill: `_create_array(_reconstruct, (ndarray, (0,), b'b'),
    /// (1, shape, dtype, False, bytes), None)`, up to the dtype.
    fn array_start(&mut self, shape: &[i64]) {
        self.get(CREATE_ARRAY);
        self.check();
        self.op(MARK);
        self.get(NDARRAY);
        self.get(B);
        self.int(0);
        self.op(TUPLE1);
        self.memoize();
        self.get(SIMPLE_B);
        self.op(TUPLE3);
        self.memoize();
        self.check();
        self.op(MARK);
        self.int(1);
        for dimension in shape {
            self.int(*dimension);
        }
        self.op(match shape.len() {
            1 => TUPLE1,
            2 => TUPLE2,
            _ => TUPLE3,
        });
        self.memoize();
    }

    fn array_end(&mut self, data: &[u8]) {
        self.falsity();
        self.bytes(data);
        self.op(TUPLE);
        self.memoize();
        self.none();
        self.op(TUPLE);
        self.memoize();
        self.op(REDUCE);
        self.memoize();
    }

    fn array(&mut self, shape: &[i64], dtype: u8, data: &[u8]) {
        self.array_start(shape);
        self.get(dtype);
        self.array_end(data);
    }

    /// `_batch_appends`, with the list's elements written by `item`.
    fn appends(&mut self, count: usize, mut item: impl FnMut(&mut Writer, usize)) {
        let mut at = 0;
        loop {
            let batch = (count - at).min(1000);
            if batch > 1 {
                self.op(MARK);
                for index in at..at + batch {
                    item(self, index);
                }
                self.op(APPENDS);
            } else if batch == 1 {
                item(self, at);
                self.op(APPEND);
            }
            at += batch;
            if batch < 1000 {
                return;
            }
        }
    }

    /// `save_tuple`, with the elements written by `item`.
    fn tuple(&mut self, count: usize, mut item: impl FnMut(&mut Writer, usize)) {
        if count <= 3 {
            for index in 0..count {
                item(self, index);
            }
            self.op([EMPTY_TUPLE, TUPLE1, TUPLE2, TUPLE3][count]);
        } else {
            self.check();
            self.op(MARK);
            for index in 0..count {
                item(self, index);
            }
            self.op(TUPLE);
        }
        self.memoize();
    }
}

fn little_doubles(values: impl Iterator<Item = f64>) -> Vec<u8> {
    values.flat_map(f64::to_le_bytes).collect()
}

fn little_ints(values: impl Iterator<Item = i64>) -> Vec<u8> {
    values.flat_map(i64::to_le_bytes).collect()
}

/// `Tree.nodes` as the struct array holds it: `missing_go_to_left` and its padding are zero.
fn node_bytes(nodes: &[Node]) -> Vec<u8> {
    let mut out = Vec::with_capacity(nodes.len() * 64);
    for node in nodes {
        out.extend_from_slice(&node.left_child.to_le_bytes());
        out.extend_from_slice(&node.right_child.to_le_bytes());
        out.extend_from_slice(&node.feature.to_le_bytes());
        out.extend_from_slice(&node.threshold.to_le_bytes());
        out.extend_from_slice(&node.impurity.to_le_bytes());
        out.extend_from_slice(&node.n_node_samples.to_le_bytes());
        out.extend_from_slice(&node.weighted_n_node_samples.to_le_bytes());
        out.extend_from_slice(&[0u8; 8]);
    }
    out
}

/// What follows the first tree's dictionary in each of the others.
fn later_tree(w: &mut Writer, forest: &Forest, index: usize) {
    let tree = &forest.trees[index];
    let columns = forest.statistics.len() as i64;
    w.get(REGRESSOR);
    w.check();
    w.op(EMPTY_TUPLE);
    w.op(NEWOBJ);
    w.memoize();
    w.check();
    w.op(EMPTY_DICT);
    w.memoize();
    w.op(MARK);
    for key in [128, 129, 130, 131] {
        w.get(key);
    }
    w.get(DEPTH_KEY);
    w.int(forest.max_depth);
    w.get(133);
    w.int(2);
    w.get(134);
    w.int(1);
    w.get(135);
    w.float(0.0);
    w.get(136);
    w.int(1);
    w.get(137);
    w.none();
    w.get(138);
    w.int(i64::from(tree.random_state));
    w.get(139);
    w.float(0.0);
    w.get(140);
    w.none();
    w.get(141);
    w.float(0.0);
    w.get(87);
    w.int(columns);
    w.get(164);
    w.int(1);
    w.get(165);
    w.int(1);
    w.get(166);
    w.get(TREE_CLASS);
    w.int(columns);
    w.array(&[1], I8, &1i64.to_le_bytes());
    w.int(1);
    w.op(TUPLE3);
    w.memoize();
    w.op(REDUCE);
    w.memoize();
    w.check();
    w.op(EMPTY_DICT);
    w.memoize();
    w.op(MARK);
    w.get(DEPTH_KEY);
    w.int(tree.max_depth);
    w.get(COUNT_KEY);
    let count = tree.nodes.len() as i64;
    w.int(count);
    w.get(NODES_KEY);
    w.array(&[count], NODE_DTYPE, &node_bytes(&tree.nodes));
    w.get(VALUES_KEY);
    w.array(
        &[count, 1, 1],
        F8,
        &little_doubles(tree.nodes.iter().map(|node| node.value)),
    );
    w.op(SETITEMS);
    w.op(BUILD);
    w.get(VERSION_KEY);
    w.get(VERSION);
    w.op(SETITEMS);
    w.op(BUILD);
}

/// The first tree, which is where the memo learns the shape of the others.
fn first_tree(w: &mut Writer, forest: &Forest) {
    let tree = &forest.trees[0];
    let columns = forest.statistics.len() as i64;
    let count = tree.nodes.len() as i64;
    w.raw(template::TREE_A);
    w.int(i64::from(tree.random_state));
    w.raw(template::TREE_B);
    w.int(columns);
    w.raw(template::TREE_C);
    w.int(columns);
    w.raw(template::TREE_D);
    w.int(tree.max_depth);
    w.raw(template::TREE_E);
    w.int(count);
    w.raw(template::TREE_F);
    w.int(count);
    w.raw(template::TREE_G);
    w.bytes(&node_bytes(&tree.nodes));
    w.raw(&template::TREE_H[1..]);
    w.int(count);
    w.raw(template::TREE_I);
    w.bytes(&little_doubles(tree.nodes.iter().map(|node| node.value)));
    w.raw(&template::TREE_J[1..]);
}

/// `dill.dump(scorer_lambda, f)`.
pub fn write(forest: &Forest) -> Vec<u8> {
    let columns = forest.statistics.len() as i64;
    let mut w = Writer {
        out: vec![PROTO, 4],
        frame: Vec::new(),
    };
    let random_state = i64::from(forest.hyperparameters.random_state);
    w.raw(template::CODE_START);
    w.string(SCRIPT_NAME);
    w.raw(&template::AFTER_NAME[1..]);
    w.int(columns);
    w.raw(template::IMPUTER_A);
    w.int(columns);
    w.raw(template::IMPUTER_B);
    w.bytes(&little_doubles(forest.statistics.iter().copied()));
    w.raw(&template::FOREST_A[1..]);
    w.int(forest.max_depth);
    w.raw(template::FOREST_B);
    w.int(random_state);
    w.raw(template::FOREST_C);
    w.int(forest.trees.len() as i64);
    w.raw(template::FOREST_D);
    match forest.hyperparameters.n_jobs {
        Some(jobs) => w.int(jobs),
        None => w.none(),
    }
    w.raw(template::FOREST_E);
    w.int(random_state);
    w.raw(template::FOREST_F);
    w.int(columns);
    w.raw(template::FOREST_G);
    w.int(forest.max_samples as i64);
    w.raw(template::FOREST_H);
    w.int(forest.n_samples as i64);
    w.raw(template::FOREST_I);
    w.int(forest.max_samples as i64);
    w.raw(template::FOREST_J);
    w.int(columns);
    w.raw(template::FOREST_K);
    w.appends(forest.trees.len(), |w, index| {
        if index == 0 {
            first_tree(w, forest);
        } else {
            later_tree(w, forest, index);
        }
    });

    w.string("estimators_features_");
    w.check();
    w.op(EMPTY_LIST);
    w.memoize();
    let features = little_ints(0..columns);
    w.appends(forest.trees.len(), |w, _| {
        w.array(&[columns], I8, &features)
    });
    w.string("_seeds");
    w.array(
        &[forest.seeds.len() as i64],
        I8,
        &little_ints(forest.seeds.iter().copied()),
    );
    w.string("_average_path_length_per_tree");
    w.tuple(forest.trees.len(), |w, index| {
        let lengths = &forest.average_path_lengths[index];
        w.array(
            &[lengths.len() as i64],
            F8,
            &little_doubles(lengths.iter().copied()),
        );
    });
    w.string("_decision_path_lengths");
    w.tuple(forest.trees.len(), |w, index| {
        let depths = &forest.decision_path_lengths[index];
        w.array(
            &[depths.len() as i64],
            I8,
            &little_ints(depths.iter().copied()),
        );
    });
    w.string("offset_");
    w.float(-0.5);
    w.get(VERSION_KEY);
    w.get(VERSION);
    w.op(SETITEMS);
    w.op(BUILD);
    w.op(TUPLE3);
    w.memoize();
    w.op(REDUCE);
    w.op(POP);

    // The third cell: `annotation_names_i`.
    w.get(42);
    w.get(62);
    w.get(43);
    let longest = forest
        .annotation_names
        .iter()
        .map(|name| name.chars().count())
        .max()
        .unwrap_or(0);
    w.array_start(&[forest.annotation_names.len() as i64]);
    w.get(DTYPE);
    w.string(&format!("U{longest}"));
    w.falsity();
    w.truth();
    w.op(TUPLE3);
    w.memoize();
    w.op(REDUCE);
    w.memoize();
    w.check();
    w.op(MARK);
    w.int(3);
    w.get(LITTLE);
    w.none();
    w.none();
    w.none();
    w.int(longest as i64 * 4);
    w.int(4);
    w.int(8);
    w.op(TUPLE);
    w.memoize();
    w.op(BUILD);
    let mut names = Vec::with_capacity(forest.annotation_names.len() * longest * 4);
    for name in &forest.annotation_names {
        let mut written = 0;
        for character in name.chars() {
            names.extend_from_slice(&u32::from(character).to_le_bytes());
            written += 1;
        }
        names.resize(names.len() + (longest - written) * 4, 0);
    }
    w.array_end(&names);
    w.op(TUPLE3);
    w.memoize();
    w.op(REDUCE);
    w.op(POP);
    w.op(TUPLE3);
    w.memoize();
    w.op(REDUCE);
    w.op(POP);
    w.op(STOP);
    w.commit(true);
    w.out
}

/// A value on the interpreter's stack. Nothing it names is ever looked up.
#[derive(Debug, Clone)]
enum Value {
    Mark,
    None,
    Bool,
    Int(i64),
    Float,
    Str(Rc<String>),
    Bytes(Rc<Vec<u8>>),
    Tuple(Rc<Vec<Value>>),
    List(Rc<RefCell<Vec<Value>>>),
    Dict(Rc<RefCell<Vec<(Value, Value)>>>),
    Global(Rc<String>),
    Call(Rc<Call>),
}

#[derive(Debug)]
struct Call {
    callable: Value,
    arguments: Value,
    state: RefCell<Option<Value>>,
}

impl Value {
    fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Int(value) => Some(*value),
            _ => None,
        }
    }

    fn call_of(&self, name: &str) -> Option<&Call> {
        match self {
            Value::Call(call) => match &call.callable {
                Value::Global(global) if global.as_str() == name => Some(call),
                _ => None,
            },
            _ => None,
        }
    }

    /// The value of a dictionary's string key.
    fn get(&self, key: &str) -> Option<Value> {
        let Value::Dict(entries) = self else {
            return None;
        };
        entries.borrow().iter().find_map(|(k, v)| match k {
            Value::Str(text) if text.as_str() == key => Some(v.clone()),
            _ => None,
        })
    }

    fn elements(&self) -> Option<Vec<Value>> {
        match self {
            Value::Tuple(values) => Some(values.as_ref().clone()),
            Value::List(values) => Some(values.borrow().clone()),
            _ => None,
        }
    }

    /// An instance's dictionary, from `NEWOBJ` and `BUILD`, or a `Tree`'s from `REDUCE`.
    fn state(&self) -> Option<Value> {
        match self {
            Value::Call(call) => call.state.borrow().clone(),
            _ => None,
        }
    }

    /// `_create_array(...)`'s shape and data.
    fn array(&self) -> Option<(Vec<i64>, Rc<Vec<u8>>, Value)> {
        let call = self.call_of("dill._dill._create_array")?;
        let arguments = call.arguments.elements()?;
        let state = arguments.get(2)?.elements()?;
        let shape = state
            .get(1)?
            .elements()?
            .iter()
            .map(Value::as_i64)
            .collect::<Option<Vec<i64>>>()?;
        let Value::Bytes(data) = state.get(4)? else {
            return None;
        };
        Some((shape, data.clone(), state.get(2)?.clone()))
    }
}

/// Interprets a pickle, keeping every array and instance that was made along the way.
fn interpret(bytes: &[u8]) -> Option<Vec<Value>> {
    let mut stack: Vec<Value> = Vec::new();
    let mut memo: HashMap<u32, Value> = HashMap::new();
    let mut made = Vec::new();
    let mut at = 0usize;
    fn take<'a>(bytes: &'a [u8], at: &mut usize, n: usize) -> Option<&'a [u8]> {
        let slice = bytes.get(*at..at.checked_add(n)?)?;
        *at += n;
        Some(slice)
    }
    fn until_mark(stack: &mut Vec<Value>) -> Option<Vec<Value>> {
        let mark = stack.iter().rposition(|v| matches!(v, Value::Mark))?;
        let items = stack.split_off(mark + 1);
        stack.pop();
        Some(items)
    }
    loop {
        let op = *take(bytes, &mut at, 1)?.first()?;
        match op {
            PROTO => {
                take(bytes, &mut at, 1)?;
            }
            FRAME => {
                take(bytes, &mut at, 8)?;
            }
            MARK => stack.push(Value::Mark),
            NONE => stack.push(Value::None),
            NEWTRUE | NEWFALSE => stack.push(Value::Bool),
            BININT1 => stack.push(Value::Int(i64::from(take(bytes, &mut at, 1)?[0]))),
            BININT2 => {
                let b = take(bytes, &mut at, 2)?;
                stack.push(Value::Int(i64::from(u16::from_le_bytes([b[0], b[1]]))));
            }
            BININT => {
                let b = take(bytes, &mut at, 4)?;
                stack.push(Value::Int(i64::from(i32::from_le_bytes(
                    b.try_into().ok()?,
                ))));
            }
            LONG1 => {
                let n = usize::from(take(bytes, &mut at, 1)?[0]);
                let b = take(bytes, &mut at, n)?;
                if n > 8 {
                    return None;
                }
                let fill = if b.last().is_some_and(|top| top & 0x80 != 0) {
                    0xff
                } else {
                    0
                };
                let mut wide = [fill; 8];
                wide[..n].copy_from_slice(b);
                stack.push(Value::Int(i64::from_le_bytes(wide)));
            }
            BINFLOAT => {
                take(bytes, &mut at, 8)?;
                stack.push(Value::Float);
            }
            SHORT_BINUNICODE => {
                let n = usize::from(take(bytes, &mut at, 1)?[0]);
                let b = take(bytes, &mut at, n)?;
                stack.push(Value::Str(Rc::new(String::from_utf8(b.to_vec()).ok()?)));
            }
            BINUNICODE => {
                let b = take(bytes, &mut at, 4)?;
                let n = u32::from_le_bytes(b.try_into().ok()?) as usize;
                let b = take(bytes, &mut at, n)?;
                stack.push(Value::Str(Rc::new(String::from_utf8(b.to_vec()).ok()?)));
            }
            SHORT_BINBYTES => {
                let n = usize::from(take(bytes, &mut at, 1)?[0]);
                stack.push(Value::Bytes(Rc::new(take(bytes, &mut at, n)?.to_vec())));
            }
            BINBYTES => {
                let b = take(bytes, &mut at, 4)?;
                let n = u32::from_le_bytes(b.try_into().ok()?) as usize;
                stack.push(Value::Bytes(Rc::new(take(bytes, &mut at, n)?.to_vec())));
            }
            EMPTY_TUPLE => stack.push(Value::Tuple(Rc::new(Vec::new()))),
            EMPTY_LIST => stack.push(Value::List(Rc::new(RefCell::new(Vec::new())))),
            EMPTY_DICT => stack.push(Value::Dict(Rc::new(RefCell::new(Vec::new())))),
            TUPLE => {
                let items = until_mark(&mut stack)?;
                stack.push(Value::Tuple(Rc::new(items)));
            }
            TUPLE1 | TUPLE2 | TUPLE3 => {
                let n = usize::from(op - TUPLE1) + 1;
                let items = stack.split_off(stack.len().checked_sub(n)?);
                stack.push(Value::Tuple(Rc::new(items)));
            }
            APPEND => {
                let item = stack.pop()?;
                let Value::List(list) = stack.last()? else {
                    return None;
                };
                list.borrow_mut().push(item);
            }
            APPENDS => {
                let items = until_mark(&mut stack)?;
                let Value::List(list) = stack.last()? else {
                    return None;
                };
                list.borrow_mut().extend(items);
            }
            SETITEM => {
                let value = stack.pop()?;
                let key = stack.pop()?;
                let Value::Dict(dict) = stack.last()? else {
                    return None;
                };
                dict.borrow_mut().push((key, value));
            }
            SETITEMS => {
                let items = until_mark(&mut stack)?;
                if items.len() % 2 != 0 {
                    return None;
                }
                let Value::Dict(dict) = stack.last()? else {
                    return None;
                };
                let mut dict = dict.borrow_mut();
                for pair in items.chunks(2) {
                    dict.push((pair[0].clone(), pair[1].clone()));
                }
            }
            MEMOIZE => {
                memo.insert(memo.len() as u32, stack.last()?.clone());
            }
            BINGET => {
                let index = u32::from(take(bytes, &mut at, 1)?[0]);
                stack.push(memo.get(&index)?.clone());
            }
            LONG_BINGET => {
                let b = take(bytes, &mut at, 4)?;
                let index = u32::from_le_bytes(b.try_into().ok()?);
                stack.push(memo.get(&index)?.clone());
            }
            GLOBAL => {
                let mut read_line = || -> Option<String> {
                    let end = bytes.get(at..)?.iter().position(|b| *b == b'\n')?;
                    let line = String::from_utf8(bytes[at..at + end].to_vec()).ok()?;
                    at += end + 1;
                    Some(line)
                };
                let module = read_line()?;
                let name = read_line()?;
                stack.push(Value::Global(Rc::new(format!("{module}.{name}"))));
            }
            STACK_GLOBAL => {
                let Value::Str(name) = stack.pop()? else {
                    return None;
                };
                let Value::Str(module) = stack.pop()? else {
                    return None;
                };
                stack.push(Value::Global(Rc::new(format!("{module}.{name}"))));
            }
            REDUCE | NEWOBJ => {
                let arguments = stack.pop()?;
                let callable = stack.pop()?;
                let value = Value::Call(Rc::new(Call {
                    callable,
                    arguments,
                    state: RefCell::new(None),
                }));
                made.push(value.clone());
                stack.push(value);
            }
            BUILD => {
                let state = stack.pop()?;
                if let Value::Call(call) = stack.last()? {
                    *call.state.borrow_mut() = Some(state);
                }
            }
            POP => {
                stack.pop()?;
            }
            STOP => return Some(made),
            _ => return None,
        }
    }
}

fn doubles_of(data: &[u8]) -> Vec<f64> {
    data.chunks_exact(8)
        .map(|b| f64::from_le_bytes(b.try_into().expect("eight bytes")))
        .collect()
}

fn ints_of(data: &[u8]) -> Vec<i64> {
    data.chunks_exact(8)
        .map(|b| i64::from_le_bytes(b.try_into().expect("eight bytes")))
        .collect()
}

/// A numpy `U` array's strings: UTF-32 at the dtype's width, padded with zeros.
fn strings_of(data: &[u8], width: usize) -> Option<Vec<String>> {
    if width == 0 {
        return Some(Vec::new());
    }
    data.chunks_exact(width * 4)
        .map(|element| {
            element
                .chunks_exact(4)
                .map(|b| u32::from_le_bytes(b.try_into().expect("four bytes")))
                .take_while(|code| *code != 0)
                .map(char::from_u32)
                .collect()
        })
        .collect()
}

/// A scorer the default script wrote at scikit-learn 1.3.2, or `None` for anything else.
pub fn read(bytes: &[u8]) -> Option<Scorer> {
    let made = interpret(bytes)?;
    let state_of = |class: &str| -> Option<Value> {
        made.iter()
            .find(|value| value.call_of(class).is_some())?
            .state()
    };
    let imputer = state_of("sklearn.impute._base.SimpleImputer")?;
    let forest = state_of("sklearn.ensemble._iforest.IsolationForest")?;
    match forest.get("_sklearn_version")? {
        Value::Str(version) if version.as_str() == "1.3.2" => {}
        _ => return None,
    }
    let (_, statistics, _) = imputer.get("statistics_")?.array()?;
    let mut trees = Vec::new();
    for estimator in forest.get("estimators_")?.elements()? {
        let state = estimator.state()?;
        let tree = state.get("tree_")?.state()?;
        let (shape, nodes, _) = tree.get("nodes")?.array()?;
        if nodes.len() != shape.first().copied()? as usize * 64 {
            return None;
        }
        let column = |offset: usize| -> Vec<[u8; 8]> {
            nodes
                .chunks_exact(64)
                .map(|node| node[offset..offset + 8].try_into().expect("eight bytes"))
                .collect()
        };
        trees.push(ScoringTree {
            left_child: column(0).into_iter().map(i64::from_le_bytes).collect(),
            right_child: column(8).into_iter().map(i64::from_le_bytes).collect(),
            feature: column(16).into_iter().map(i64::from_le_bytes).collect(),
            threshold: column(24).into_iter().map(f64::from_le_bytes).collect(),
        });
    }
    let mut average_path_lengths = Vec::new();
    for array in forest.get("_average_path_length_per_tree")?.elements()? {
        average_path_lengths.push(doubles_of(&array.array()?.1));
    }
    let mut decision_path_lengths = Vec::new();
    for array in forest.get("_decision_path_lengths")?.elements()? {
        decision_path_lengths.push(ints_of(&array.array()?.1));
    }
    if average_path_lengths.len() != trees.len() || decision_path_lengths.len() != trees.len() {
        return None;
    }
    // The last array of strings is the third cell's `annotation_names_i`.
    let mut annotation_names = None;
    for value in &made {
        if let Some((_, data, dtype)) = value.array() {
            if let Some(dtype) = dtype.call_of("numpy.dtype") {
                if let Value::Str(kind) = &dtype.arguments.elements()?[0] {
                    if let Some(width) = kind.strip_prefix('U') {
                        annotation_names = strings_of(&data, width.parse().ok()?);
                    }
                }
            }
        }
    }
    let max_samples = forest.get("_max_samples")?.as_i64()?;
    if !(1..=256).contains(&max_samples) {
        return None;
    }
    Some(Scorer {
        annotation_names: annotation_names?,
        statistics: doubles_of(&statistics),
        max_samples,
        trees,
        average_path_lengths,
        decision_path_lengths,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isolation_forest::{fit, score_samples, Hyperparameters};

    fn forest(trees: usize) -> Forest {
        let rows: Vec<Vec<f64>> = (0..40)
            .map(|i| vec![f64::from(i), f64::from(i * 7 % 13)])
            .collect();
        fit(
            &["AA".to_string(), "B".to_string()],
            &rows,
            &Hyperparameters {
                random_state: 0,
                n_estimators: trees,
                n_jobs: None,
            },
        )
        .unwrap()
    }

    #[test]
    fn a_written_scorer_reads_back_and_scores_alike() {
        for trees in [1, 2, 3, 4, 100] {
            let forest = forest(trees);
            let read_back = read(&write(&forest)).expect("a scorer");
            assert_eq!(read_back, forest.scorer(), "{trees} trees");
        }
        let names = ["AA".to_string(), "B".to_string()];
        let rows = vec![vec![1.0, 2.0], vec![f64::NAN, 30.0]];
        let scorer = forest(10).scorer();
        assert_eq!(
            score_samples(&scorer, &names, &rows),
            score_samples(&read(&write(&forest(10))).unwrap(), &names, &rows)
        );
    }

    #[test]
    fn what_is_not_this_scripts_scorer_is_not_read() {
        assert!(read(b"no pickle\n").is_none());
        assert!(read(&[]).is_none());
        // A valid pickle of nothing the port knows.
        assert!(read(&[PROTO, 4, NONE, STOP]).is_none());
    }

    #[test]
    fn a_frame_closes_at_the_first_save_past_the_target() {
        let forest = forest(100);
        let bytes = write(&forest);
        assert_eq!(bytes[2], FRAME);
    }
}
