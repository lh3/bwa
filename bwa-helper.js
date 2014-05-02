/*****************************************************************
 * The K8 Javascript interpreter is required to run this script. *
 *                                                               *
 *   Source code: https://github.com/attractivechaos/k8          *
 *   Binary: https://sourceforge.net/projects/lh3/files/         *
 *****************************************************************/

/******************
 *** From k8.js ***
 ******************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/************************
 *** command markovlp ***
 ************************/

function bwa_markOvlp(args)
{
	var c, min_aln_ratio = .9, min_ext = 50;
	while ((c = getopt(args, "r:e:")) != null) {
		if (c == 'r') min_aln_ratio = parseFloat(getopt.arg);
		else if (c == 'e') min_ext = parseInt(getopt.arg);
	}

	var file = args.length == getopt.ind? new File() : new File(args[getopt.ind]);
	var buf = new Bytes();
	var dir4 = ['>>', '><', '<>', '<<'];

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t")
			for (var i = 0; i < t.length; ++i)
				if (i != 0 && i != 4)
					t[i] = parseInt(t[i]);
		var el, a1, a2, e1, e2, r1, r2; // a: aligned length; e: extended length; r: remaining length
		e2 = a2 = t[7] - t[6];
		if (t[2] < t[3]) { // forward-forward match
			e1 = a1 = t[3] - t[2];
			r1 = t[2] - t[6];
			r2 = (t[5] - t[7]) - (t[1] - t[3]);
			el = r1 > 0? t[6] : t[2];
			el += r2 > 0? t[1] - t[3] : t[5] - t[7];
		} else { // reverse-forward match
			e1 = a1 = t[2] - t[3];
			r1 = (t[1] - t[2]) - t[6];
			r2 = (t[5] - t[7]) - t[3];
			el = r1 > 0? t[6] : t[1] - t[2];
			el += r2 > 0? t[3] : t[5] - t[7];
		}
		e1 += el; e2 += el;
		var type;
		if (a1 / e1 >= min_aln_ratio && a2 / e2 >= min_aln_ratio) {
			if ((r1 >= min_ext && r2 >= min_ext) || (r1 <= -min_ext && r2 <= -min_ext)) { // suffix-prefix match
				var d = t[2] < t[3]? 0 : 2;
				if (r1 < 0) d ^= 3; // reverse the direction
				type = 'O' + dir4[d];
			} else type = 'C' + (e1 / t[1] > e2 / t[5]? 1 : 2);
		} else type = 'I'; // internal local match; not a suffix-prefix match
		//print(t[1], e1, a1, t[5], e2, a2);
		print(type, buf);
	}

	buf.destroy();
	file.close();
}

/***********************
 *** command pas2bed ***
 ***********************/

function bwa_pas2reg(args)
{
	var file = args.length? new File(args[0]) : new File();
	var buf = new Bytes();

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0] == t[4]) continue;
		if (parseInt(t[2]) < parseInt(t[3])) print(t[0], t[1], t[2], t[3], t[8]);
		else print(t[0], t[1], t[3], t[2], t[8]);
		print(t[4], t[5], t[6], t[7], t[8]);
	}

	buf.destroy();
	file.close();
}

/*******************
 * command sam2pas *
 *******************/

function bwa_sam2pas(args)
{
	var file = args.length == 0? new File() : new File(args[0]);
	var buf = new Bytes();
	var seq_dict = {};

	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var m;
		if (/^@SQ/.test(line)) {
			var name = null, len = null;
			if ((m = /\tSN:(\S+)/.exec(line)) != null) name = m[1];
			if ((m = /\tLN:(\S+)/.exec(line)) != null) len = parseInt(m[1]);
			if (name != null && len != null) seq_dict[name] = len;
		}
		if (/^@/.test(line)) continue;
		var t = line.split("\t");
		var pos = parseInt(t[3]) - 1;
		var x = 0, y = 0, i = 0, clip = [0, 0], n_ins = 0, n_del = 0, o_ins = 0, o_del = 0, n_M = 0;
		var re = /(\d+)([MIDSH])/g;
		while ((m = re.exec(t[5])) != null) {
			var l = parseInt(m[1]);
			if (m[2] == 'M') x += l, y += l, n_M += l;
			else if (m[2] == 'I') y += l, n_ins += l, ++o_ins;
			else if (m[2] == 'D') x += l, n_del += l, ++o_del;
			else if (m[2] == 'S' || m[2] == 'H')
				clip[i == 0? 0 : 1] = l;
			++i;
		}
		var is_rev = (parseInt(t[1]) & 16)? true : false;
		var misc = 'mapQ=' + t[4] + ';';
		var usc = 1;
		if ((m = /\tNM:i:(\d+)/.exec(line)) != null) {
			var NM = parseInt(m[1]);
			var diff = (NM / (n_M + n_ins + n_del)).toFixed(3);
			misc += 'diff=' + diff + ';n_mis=' + (NM - n_del - n_ins) + ';';
		}
		misc += 'n_del='+n_del+';n_ins='+n_ins+';o_del='+o_del+';o_ins='+o_ins + ';';
		if ((m = /\tAS:i:(\d+)/.exec(line)) != null) {
			misc += 'AS='+m[1] + ';';
			usc = (parseInt(m[1]) / (x > y? x : y)).toFixed(3);
		}
		if ((m = /\tXS:i:(\d+)/.exec(line)) != null) misc += 'XS='+m[1] + ';';
		var len = y + clip[0] + clip[1];
		var z = [t[0], len, clip[0], clip[0] + y, t[2], seq_dict[t[2]], pos, pos + x, usc, misc];
		if (parseInt(t[1]) & 16) z[2] = clip[1] + y, z[3] = clip[1];
		print(z.join("\t"));
	}

	buf.destroy();
	file.close();
}

/***********************
 *** command reg2cut ***
 ***********************/

function bwa_reg2cut(args)
{
	var c, min_usc = 0.5, min_ext = 100, min_len = 5000, cut = 250;
	while ((c = getopt(args, "s:e:l:c:")) != null) {
		if (c == 's') min_usc = parseFloat(getopt.arg);
		else if (c == 'e') min_ext = parseInt(getopt.arg);
		else if (c == 'l') min_len = parseInt(getopt.arg);
		else if (c == 'c') cut = parseInt(getopt.arg);
	}

	var file = args.length == getopt.ind? new File () : new File(args[getopt.ind]);
	var buf = new Bytes();

	function print_bed() {
		for (var i = 0; i < a.length; ++i) {
			var start = a[i][0] - cut > 0? a[i][0] : 0;
			var end   = a[i][1] + cut < last_len? a[i][1] : last_len;
			if (end - start >= min_len) print(last_chr, start, end);
		}
	}

	var last_chr = null, last_len = null, max_c_usc = 0, start = 0, end = 0;
	var a = [];
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		t[1] = parseInt(t[1]);
		t[2] = parseInt(t[2]);
		t[3] = parseInt(t[3]);
		t[4] = parseFloat(t[4]);
		var is_contained = t[2] < min_ext && t[1] - t[3] < min_ext? true : false;
		if (t[3] - t[2] < cut<<1) continue;
		t[2] += cut; t[3] -= cut;
		if (t[0] != last_chr) {
			a.push([start, end]);
			if (last_chr != null && max_c_usc < min_usc) print_bed();
			last_chr = t[0]; last_len = t[1]; start = t[2]; end = t[3];
			max_c_usc = is_contained? t[4] : 0;
			a.length = 0;
		} else {
			if (is_contained)
				max_c_usc = max_c_usc > t[4]? max_c_usc : t[4];
			if (t[4] < min_usc) continue;
			if (t[2] > end) {
				a.push([start, end]);
				start = t[2];
				end = end > t[3]? end : t[3];
			} else end = end > t[3]? end : t[3];
		}
	}
	a.push([start, end]);
	if (max_c_usc < min_usc) print_bed(); // the last sequence

	buf.destroy();
	file.close();
}

function bwa_shortname(args)
{
	var file = args.length? new File(args[0]) : new File();
	var buf = new Bytes();

	var re = /(\S+)\/(\d+)_(\d+)((:\d+-\d+)+)/g;
	var re2 = /:(\d+)-(\d+)/g;
	while (file.readline(buf) >= 0) {
		var match, match2;
		var line = buf.toString();
		var x = [];
		while ((match = re.exec(line)) != null) {
			var start = parseInt(match[2]), len = parseInt(match[3]) - start;
			while ((match2 = re2.exec(match[4])) != null) {
				var a = parseInt(match2[1]) - 1;
				var b = parseInt(match2[2]);
				start += a; len = b - a;
			}
			x.push([match[0], match[1] + '/' + start.toString() + '_' + (start+len).toString()]);
		}
		for (var i = 0; i < x.length; ++i)
			line = line.replace(x[i][0], x[i][1]);
		print(line);
	}

	buf.destroy();
	file.close();
}

/*******************
 * Command gff2sam *
 *******************/

function bwa_gff2sam(args)
{
	if (args.length < 2) {
		print("Usage: k8 bwa-helper.js <aln.gff> <query-length.txt>");
		exit(1);
	}

	var file = new File(args[1]);
	var buf = new Bytes();
	var len = {};

	while (file.readline(buf) >= 0) {
		var t = buf.toString().split(/\s+/);
		len[t[0]] = parseInt(t[1]);
	}
	file.close();

	file = new File(args[0]);
	var re_cigar = /([MID])(\d+)/g;
	var lineno = 0;
	while (file.readline(buf) >= 0) {
		++lineno;
		var t = buf.toString().split("\t");
		var m = /Target=(\S+)\s+(\d+)\s+(\d+)\s+([+-])/.exec(t[8]);
		if (m == null) {
			warn("WARNING: skipped line "+lineno+" due to the lack of Target.");
			continue;
		}
		var qname = m[1];
		var flag = t[6] == m[4]? 0 : 16;
		var qb = parseInt(m[2]) - 1, qe = parseInt(m[3]), qlen = len[qname];
		if (qlen == null)
			throw Error("Sequence "+qname+" is not present in the query-length.txt");
		var clip5 = qb, clip3 = qlen - qe;
		if (flag&16) clip5 ^= clip3, clip3 ^= clip5, clip5 ^= clip3; // swap

		m = /Gap\s*=\s*(([MID]\d+\s*)+)/.exec(t[8]);
		var cigar = clip5? clip5 + 'S' : '';
		var n_ins = 0, n_del = 0, n_match = 0, NM = null;
		if (m) {
			var mc;
			while ((mc = re_cigar.exec(m[1])) != null) {
				var l = parseInt(mc[2]);
				cigar += mc[2] + mc[1];
				if (mc[1] == 'I') n_ins += l;
				else if (mc[1] == 'D') n_del += l;
				else if (mc[1] == 'M') n_match += l;
			}
			if (n_ins + n_match != qe - qb || n_del + n_match != parseInt(t[4]) - parseInt(t[3]) + 1)
				throw Error("Inconsistent CIGAR at line "+lineno);
		} else { // ungapped alignment
			var tb = parseInt(t[3]) - 1, te = parseInt(t[4]);
			if (te - tb != qe - qb) {
				warn("WARNING: line "+lineno+" should contain gaps, but lacks Gap. Skipped.\n");
			} else cigar = (qe - qb) + 'M';
		}
		if (clip3) cigar += clip3 + 'S';
		if ((m = /num_mismatch=(\d+)/.exec(t[8])) != null)
			NM = parseInt(m[1]) + n_ins + n_del;
		var out = [qname, flag, t[0], t[3], 255, cigar, '*', 0, 0, '*', '*'];
		if (NM != null) out.push('NM:i:' + NM);
		print(out.join("\t"));
	}
	buf.destroy();
	file.close();
}

/*********************
 *** Main function ***
 *********************/

function main(args)
{
	if (args.length == 0) {
		print("\nUsage:    k8 bwa-helper.js <command> [arguments]\n");
		print("Commands: sam2pas      convert SAM to pairwise alignment summary format (PAS)");
		print("          pas2reg      extract covered regions");
		print("          reg2cut      regions to extract for the 2nd round bwa-mem");
		print("          markovlp     identify bi-directional overlaps");
		print("          gff2sam      convert GFF3 alignment to SAM");
		print("          shortname    shorten sequence name after subseq (PacBio read names only)");
		print("");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'sam2pas') bwa_sam2pas(args);
	else if (cmd == 'gff2sam') bwa_gff2sam(args);
	else if (cmd == 'markovlp') bwa_markOvlp(args);
	else if (cmd == 'pas2reg') bwa_pas2reg(args);
	else if (cmd == 'reg2cut') bwa_reg2cut(args);
	else if (cmd == 'shortname') bwa_shortname(args);
	else warn("Unrecognized command");
}

main(arguments);
