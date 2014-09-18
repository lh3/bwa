/*****************************************************************
 * The K8 Javascript interpreter is required to run this script. *
 *                                                               *
 *   Source code: https://github.com/attractivechaos/k8          *
 *   Binary: http://sourceforge.net/projects/lh3/files/k8/       *
 *                                                               *
 * Data file used for generating GRCh38 ALT alignments:          *
 *                                                               *
 *   http://sourceforge.net/projects/bio-bwa/files/              *
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

function obj2str(o)
{
	if (typeof(o) != 'object') {
		return o.toString();
	} else if (o == null) {
		return "null";
	} else if (Array.isArray(o)) {
		var s = "[";
		for (var i = 0; i < o.length; ++i) {
			if (i) s += ',';
			s += obj2str(o[i]);
		}
		return s + "]";
	} else {
		var i = 0, s = "{";
		for (var key in o) {
			if (i++) s += ',';
			s += key + ":";
			s += obj2str(o[key]);
		}
		return s + "}";
	}
}

Bytes.prototype.reverse = function()
{
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[i];
		this[i] = this[this.length - i - 1];
		this[this.length - i - 1] = tmp;
	}
}

Bytes.prototype.revcomp = function()
{
	if (Bytes.rctab == null) {
		var s1 = 'WSATUGCYRKMBDHVNwsatugcyrkmbdhvn';
		var s2 = 'WSTAACGRYMKVHDBNwstaacgrymkvhdbn';
		Bytes.rctab = [];
		for (var i = 0; i < 256; ++i) Bytes.rctab[i] = 0;
		for (var i = 0; i < s1.length; ++i)
			Bytes.rctab[s1.charCodeAt(i)] = s2.charCodeAt(i);
	}
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[this.length - i - 1];
		this[this.length - i - 1] = Bytes.rctab[this[i]];
		this[i] = Bytes.rctab[tmp];
	}
	if (this.length>>1)
		this[this.length>>1] = Bytes.rctab[this[this.length>>1]];
}

var re_cigar = /(\d+)([MIDSHN])/g;

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


/******************************
 *** Generate ALT alignment ***
 ******************************/

function intv_ovlp(intv, bits) // interval index
{
	if (typeof bits == "undefined") bits = 13;
	intv.sort(function(a,b) {return a[0]-b[0];});
	// create the index
	var idx = [], max = 0;
	for (var i = 0; i < intv.length; ++i) {
		var b = intv[i][0]>>bits;
		var e = (intv[i][1]-1)>>bits;
		if (b != e) {
			for (var j = b; j <= e; ++j)
				if (idx[j] == null) idx[j] = i;
		} else if (idx[b] == null) idx[b] = i;
		max = max > e? max : e;
	}
	// closure
	return function(_b, _e) {
		var x = _b >> bits;
		if (x > max) return [];
		var off = idx[x];
		if (off == null) {
			var i;
			for (i = ((_e - 1) >> bits) - 1; i >= 0; --i)
				if (idx[i] != null) break;
			off = i < 0? 0 : idx[i];
		}
		var ovlp = [];
		for (var i = off; i < intv.length && intv[i][0] < _e; ++i)
			if (intv[i][1] > _b) ovlp.push(intv[i]);
		return ovlp;
	}
}

function cigar2pos(cigar, pos) // given a pos on ALT and the ALT-to-REF CIGAR, find the pos on REF
{
	var x = 0, y = 0;
	for (var i = 0; i < cigar.length; ++i) {
		var op = cigar[i][0], len = cigar[i][1];
		if (op == 'M') {
			if (y <= pos && pos < y + len)
				return x + (pos - y);
			x += len, y += len;
		} else if (op == 'D') {
			x += len;
		} else if (op == 'I') {
			if (y <= pos && pos < y + len)
				return x;
			y += len;
		} else if (op == 'S' || op == 'H') {
			if (y <= pos && pos < y + len)
				return -1;
			y += len;
		}
	}
	return -1;
}

function parse_hit(s, opt) // parse a hit. s looks something like ["chr1", "+12345", "100M", 5]
{
	var h = {};
	h.ctg = s[0];
	h.start = parseInt(s[1].substr(1)) - 1;
	h.rev = (s[1].charAt(0) == '-');
	h.cigar = s[2];
	h.NM = parseInt(s[3]);
	h.hard = false;
	var m, l_ins, n_ins, l_del, n_del, l_match, l_skip, l_clip;
	l_ins = l_del = n_ins = n_del = l_match = l_skip = l_clip = 0;
	while ((m = re_cigar.exec(h.cigar)) != null) {
		var l = parseInt(m[1]);
		if (m[2] == 'M') l_match += l;
		else if (m[2] == 'D') ++n_del, l_del += l;
		else if (m[2] == 'I') ++n_ins, l_ins += l;
		else if (m[2] == 'N') l_skip += l;
		else if (m[2] == 'H' || m[2] == 'S') {
			l_clip += l;
			if (m[2] == 'H') h.hard = true;
		}
	}
	h.end = h.start + l_match + l_del + l_skip;
	h.NM = h.NM > l_del + l_ins? h.NM : l_del + l_ins;
	h.score = Math.floor((opt.a * l_match - (opt.a + opt.b) * (h.NM - l_del - l_ins) - opt.o * (n_del + n_ins) - opt.e * (l_del + l_ins)) / opt.a + .499);
	h.l_query = l_match + l_ins + l_clip;
	return h;
}

// read the ALT-to-REF alignment and generate the index
function read_ALT_sam(fn)
{
	var intv = {};
	var file = new File(fn);
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t");
		if (line.charAt(0) == '@') continue;
		var flag = parseInt(t[1]);
		var m, cigar = [], l_qaln = 0, l_qclip = 0;
		while ((m = re_cigar.exec(t[5])) != null) {
			var l = parseInt(m[1]);
			cigar.push([m[2] != 'H'? m[2] : 'S', l]); // convert hard clip to soft clip
			if (m[2] == 'M' || m[2] == 'I') l_qaln += l;
			else if (m[2] == 'S' || m[2] == 'H') l_qclip += l;
		}
		var j = flag&16? cigar.length-1 : 0;
		var start = cigar[j][0] == 'S'? cigar[j][1] : 0;
		if (intv[t[0]] == null) intv[t[0]] = [];
		intv[t[0]].push([start, start + l_qaln, l_qaln + l_qclip, t[2], flag&16? true : false, parseInt(t[3]) - 1, cigar]);
		//print(start, start + l_qaln, t[2], flag&16? true : false, parseInt(t[3]), cigar);
	}
	buf.destroy();
	file.close();
	// create the interval index
	var idx = {};
	for (var ctg in intv)
		idx[ctg] = intv_ovlp(intv[ctg]);
	return idx;
}

function bwa_altgen(args)
{
	var c, opt = { a:1, b:4, o:6, e:1, verbose:3 };

	while ((c = getopt(args, 'v:')) != null) {
		if (c == 'v') opt.verbose = parseInt(getopt.arg);
	}

	if (args.length == getopt.ind) {
		print("Usage: k8 bwa-helper.js altgen <alt.sam> [aln.sam]");
		exit(1);
	}

	var file, buf = new Bytes();
	var aux = new Bytes(); // used for reverse and reverse complement
	var idx = read_ALT_sam(args[getopt.ind]);

	// process SAM
	file = args.length - getopt.ind >= 2? new File(args[getopt.ind+1]) : new File();
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if (line.charAt(0) == '@' || (m = /\tXA:Z:(\S+)/.exec(line)) == null) { // TODO: this does not work with PE file
			if (opt.verbose < 4) print(line);
			continue;
		}

		// parse hits
		var hits = [];
		var XA_strs = m[1].split(";");
		var NM = (m = /\tNM:i:(\d+)/.exec(line)) == null? '0' : m[1];
		var t = line.split("\t");
		var flag = parseInt(t[1]);
		var h = parse_hit([t[2], ((flag&16)?'-':'+') + t[3], t[5], NM], opt);
		if (h.hard) { // the following does not work with hard clipped SEQ
			print(line);
			continue;
		}
		hits.push(h);
		for (var i = 0; i < XA_strs.length; ++i) // hits in the XA tag
			if (XA_strs[i] != '') // as the last symbol in an XA tag is ";", the last split is an empty string
				hits.push(parse_hit(XA_strs[i].split(","), opt));

		// lift mapping positions to coordinates on the primary assembly
		var n_lifted = 0;
		for (var i = 0; i < hits.length; ++i) {
			var h = hits[i];

			if (idx[h.ctg] == null) continue;
			var a = idx[h.ctg](h.start, h.end);
			if (a == null || a.length == 0) continue;

			// find the approximate position on the primary assembly
			var lifted = [];
			for (var j = 0; j < a.length; ++j) {
				var s, e;
				if (!a[j][4]) { // ALT is mapped to the forward strand of the primary assembly
					s = cigar2pos(a[j][6], h.start);
					e = cigar2pos(a[j][6], h.end - 1) + 1;
				} else {
					s = cigar2pos(a[j][6], a[j][2] - h.end);
					e = cigar2pos(a[j][6], a[j][2] - h.start - 1) + 1;
				}
				if (s < 0 || e < 0) continue; // read is mapped to clippings in the ALT-to-chr alignment
				s += a[j][5]; e += a[j][5];
				lifted.push([a[j][3], (h.rev!=a[j][4]), s, e]);
			}
			if (lifted.length) ++n_lifted, hits[i].lifted = lifted;
		}
		if (n_lifted == 0) {
			if (opt.verbose < 4) print(line);
			continue;
		}

		// group hits
		for (var i = 0; i < hits.length; ++i) { // set keys for sorting
			if (hits[i].lifted && hits[i].lifted.length) // TODO: only the first element in lifted[] is used
				hits[i].pctg = hits[i].lifted[0][0], hits[i].pstart = hits[i].lifted[0][2], hits[i].pend = hits[i].lifted[0][3];
			else hits[i].pctg = hits[i].ctg, hits[i].pstart = hits[i].start, hits[i].pend = hits[i].end;
			hits[i].i = i; // keep the original index
		}
		hits.sort(function(a,b) { return a.pctg != b.pctg? (a.pctg < b.pctg? -1 : 1) : a.pstart - b.pstart });
		var last_chr = null, end = 0, g = -1;
		for (var i = 0; i < hits.length; ++i) {
			if (last_chr != hits[i].pctg) ++g, last_chr = hits[i].pctg, end = 0;
			else if (hits[i].pstart >= end) ++g;
			hits[i].g = g;
			end = end > hits[i].pend? end : hits[i].pend;
		}
		var reported_g = null, reported_i = null;
		for (var i = 0; i < hits.length; ++i)
			if (hits[i].i == 0)
				reported_g = hits[i].g, reported_i = i;
		var n_group0 = 0; // #hits overlapping the reported hit
		for (var i = 0; i < hits.length; ++i)
			if (hits[i].g == reported_g)
				++n_group0;
		if (n_group0 == 1) { // then keep the reported alignment and mapQ
			if (opt.verbose < 4) print(line);
			continue;
		}

		// re-estimate mapQ
		var group_max = [];
		for (var i = 0; i < hits.length; ++i) {
			var g = hits[i].g;
			if (group_max[g] == null || group_max[g][0] < hits[i].score)
				group_max[g] = [hits[i].score, g];
		}
		if (group_max.length > 1)
			group_max.sort(function(x,y) {return y[0]-x[0]});
		var mapQ;
		if (group_max[0][1] == reported_g) { // the best hit is the hit reported in SAM
			mapQ = group_max.length == 1? 60 : 6 * (group_max[0][0] - group_max[1][0]);
		} else mapQ = 0;
		mapQ = mapQ < 60? mapQ : 60;
		var ori_mapQ = parseInt(t[4]);
		mapQ = mapQ > ori_mapQ? mapQ : ori_mapQ;

		// generate lifted_str
		for (var i = 0; i < hits.length; ++i) {
			if (hits[i].lifted && hits[i].lifted.length) {
				var lifted = hits[i].lifted;
				var u = '';
				for (var j = 0; j < lifted.length; ++j)
					u += lifted[j][0] + "," + lifted[j][2] + "," + lifted[j][3] + "," + (lifted[j][1]?'-':'+') + ";";
				hits[i].lifted_str = u;
			}
		}

		// generate reversed quality and reverse-complemented sequence if necessary
		var rs = null, rq = null; // reversed quality and reverse complement sequence
		var need_rev = false;
		for (var i = 0; i < hits.length; ++i) {
			if (hits[i].g != reported_g || i == reported_i) continue;
			if (hits[i].rev != hits[reported_i].rev)
				need_rev = true;
		}
		if (need_rev) { // reverse and reverse complement
			aux.set(t[9], 0); aux.revcomp(); rs = aux.toString();
			aux.set(t[10],0); aux.reverse(); rq = aux.toString();
		}

		// print
		t[4] = mapQ;
		t.push("om:i:"+ori_mapQ);
		if (hits[reported_i].lifted_str) t.push("lt:Z:" + hits[reported_i].lifted_str);
		print(t.join("\t"));
		var cnt = 0;
		for (var i = 0; i < hits.length; ++i) {
			if (opt.verbose >= 5) print(obj2str(hits[i]));
			if (hits[i].g != reported_g || i == reported_i) continue;
			var s = [t[0], flag&0xf10, hits[i].ctg, hits[i].start+1, mapQ, hits[i].cigar, '*', 0, 0];
			// update name
			if (flag&0x40) s[0] += "/1";
			if (flag&0x80) s[0] += "/2";
			s[0] += "_" + (++cnt);
			if (hits[i].rev == hits[reported_i].rev) s.push(t[9], t[10]);
			else s.push(rs, rq);
			s.push("NM:i:" + hits[i].NM);
			if (hits[i].lifted_str) s.push("lt:Z:" + hits[i].lifted_str);
			print(s.join("\t"));
		}
	}
	file.close();

	aux.destroy();
	buf.destroy();
}

// This is in effect a simplified version of bwa_genalt().
function bwa_altlift(args)
{
	var opt = { a:1, b:4, o:6, e:1 };
	if (args.length == 0) {
		print("Usage: k8 bwa-helper.js altlift <alt-to-ref.sam> [aln.sam]");
		exit(1);
	}
	var idx = read_ALT_sam(args[0]);

	// process SAM
	var file = args.length >= 2? new File(args[1]) : new File();
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if (line.charAt(0) == '@') {
			print(line);
			continue;
		}

		var t = line.split("\t");
		var NM = (m = /\tNM:i:(\d+)/.exec(line)) == null? '0' : m[1];
		var flag = parseInt(t[1]);
		var h = parse_hit([t[2], ((flag&16)?'-':'+') + t[3], t[5], NM], opt);

		// lift mapping positions to coordinates on the primary assembly
		var a = null;
		if (idx[h.ctg] != null)
			a = idx[h.ctg](h.start, h.end);
		if (a == null) a = [];

		// find the approximate position on the primary assembly
		var lifted = [];
		for (var j = 0; j < a.length; ++j) {
			var s, e;
			if (!a[j][4]) { // ALT is mapped to the forward strand of the primary assembly
				s = cigar2pos(a[j][6], h.start);
				e = cigar2pos(a[j][6], h.end - 1) + 1;
			} else {
				s = cigar2pos(a[j][6], a[j][2] - h.end);
				e = cigar2pos(a[j][6], a[j][2] - h.start - 1) + 1;
			}
			if (s < 0 || e < 0) continue; // read is mapped to clippings in the ALT-to-chr alignment
			s += a[j][5]; e += a[j][5];
			lifted.push([a[j][3], (h.rev!=a[j][4]), s, e]);
		}
		h.lifted = lifted;

		// generate lifted_str
		if (h.lifted && h.lifted.length) {
			var lifted = h.lifted;
			var u = '';
			for (var j = 0; j < lifted.length; ++j)
				u += lifted[j][0] + "," + lifted[j][2] + "," + lifted[j][3] + "," + (lifted[j][1]?'-':'+') + ";";
			h.lifted_str = u;
		} else h.lifted_str = null;

		// print
		if (h.lifted_str) t.push("lt:Z:" + h.lifted_str);
		print(t.join("\t"));
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
		print("Commands: altlift      add lt tag to show lifted position on the primary assembly");
		print("          altgen       generate ALT alignments for ALT-unaware alignments\n");
		print("          sam2pas      convert SAM to pairwise alignment summary format (PAS)");
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
	else if (cmd == 'altgen' || cmd == 'genalt') bwa_alt(args);
	else if (cmd == 'altlift') bwa_altlift(args);
	else if (cmd == 'shortname') bwa_shortname(args);
	else warn("Unrecognized command");
}

main(arguments);
