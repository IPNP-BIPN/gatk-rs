// The fast half of the comparison: the same corpus through a clang build of the r53 sources, with
// the JNI's read loop copied literally and a fork per assembly, since an assert must not end the run.
//
//   clang -arch x86_64 -O2 -ffp-contract=off -I$R53 native_driver.c $R53/{bfc,bseq,bubble,htab,
//       kthread,ksw,mag,misc,mrope,rle,rld0,rope,unitig}.c -lz -lpthread -o driver
//   driver corpus.bin dump.txt
//
// The JNI is the oracle: this binary's memory layout is not the .so's. See README.md.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include "fml.h"

static uint32_t rd32(FILE *f, int *ok)
{
	unsigned char b[4];
	if (fread(b, 1, 4, f) != 4) { *ok = 0; return 0; }
	return (uint32_t)b[0] << 24 | b[1] << 16 | b[2] << 8 | b[3];
}

static void one(FILE *out, long idx, const char *name, int flag, int n, const char *rbuf)
{
	fml_opt_t opt;
	bseq1_t *reads = calloc(n, sizeof(bseq1_t));
	int i, j, n_utg;
	const char *p = rbuf;
	fml_utg_t *u;
	fml_opt_init(&opt);
	opt.mag_opt.flag = flag;
	for (i = 0; i < n; ++i) {
		int skip = strlen(p);
		reads[i].l_seq = skip++;
		reads[i].seq = strdup(p);
		p += skip;
		reads[i].qual = strdup(p);
		p += skip;
	}
	u = fml_assemble(&opt, n, reads, &n_utg);
	fprintf(out, "A %ld %s flag=%d reads=%d contigs=%d\n", idx, name, flag, n, n_utg);
	for (i = 0; i < n_utg; ++i) {
		fprintf(out, "C %d %d %d ", i, u[i].len, u[i].nsr);
		fwrite(u[i].seq, 1, u[i].len, out);
		fputc(' ', out);
		fwrite(u[i].cov, 1, u[i].len, out);
		fputc('\n', out);
		for (j = 0; j < u[i].n_ovlp[0] + u[i].n_ovlp[1]; ++j) {
			const fml_ovlp_t *o = &u[i].ovlp[j];
			fprintf(out, "L %d %d %d %d\n", (int)o->id, (int)o->len, !o->from, (int)o->to);
		}
	}
	fflush(out);
}

int main(int argc, char *argv[])
{
	FILE *f, *out;
	char magic[6];
	long idx = 0;
	if (argc < 3) return 2;
	// A descriptor shared with the forked child would have its offset moved when abort flushes the
	// child's streams; a memory stream has none.
	{
		FILE *disk = fopen(argv[1], "rb");
		char *all;
		long size;
		if (!disk) return 2;
		fseek(disk, 0, SEEK_END);
		size = ftell(disk);
		rewind(disk);
		all = malloc(size);
		if (!all || fread(all, 1, size, disk) != (size_t)size) return 2;
		fclose(disk);
		f = fmemopen(all, size, "rb");
	}
	out = fopen(argv[2], "wb");
	if (!f || !out || fread(magic, 1, 6, f) != 6 || memcmp(magic, "FMLC1\n", 6)) return 2;
	for (;;) {
		int ok = 1, flag, n, i, cap = 0;
		char name[256], *buf;
		size_t at = 0, room;
		uint32_t l = rd32(f, &ok);
		pid_t pid;
		int st;
		if (!ok) break;
		if (l >= sizeof(name) || fread(name, 1, l, f) != l) return 3;
		name[l] = 0;
		flag = (int)rd32(f, &ok);
		n = (int)rd32(f, &ok);
		room = 1 << 16;
		buf = calloc(room, 1);
		for (i = 0; i < 2 * n; ++i) {
			uint32_t m = rd32(f, &ok);
			if (at + m + 2 > room) {
				size_t old = room;
				while (at + m + 2 > room) room <<= 1;
				buf = realloc(buf, room);
				memset(buf + old, 0, room - old);
			}
			if (m && fread(buf + at, 1, m, f) != m) return 3;
			at += m + 1;
			cap += m + 1;
		}
		fflush(out);
		pid = fork();
		if (pid == 0) {
			if (!getenv("KEEP_STDERR")) freopen("/dev/null", "w", stderr);
			one(out, idx, name, flag, n, buf);
			_exit(0);
		}
		waitpid(pid, &st, 0);
		if (!WIFEXITED(st) || WEXITSTATUS(st) != 0)
			fprintf(out, "X %ld %s flag=%d reads=%d ABORT\n", idx, name, flag, n);
		free(buf);
		++idx;
	}
	fclose(out);
	return 0;
}
