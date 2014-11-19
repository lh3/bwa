var min_ovlp = 30;

if (arguments.length < 3) {
	print("Usage: k8 selctg.js <HLA-gene> <HLA-ALT-exons.bed> <ctg-to-ALT.sam> [min_ovlp="+min_ovlp+"]");
	exit(1);
}

if (arguments.length >= 4) min_ovlp = parseInt(arguments[3]);
var gene = arguments[0];

var buf = new Bytes();

var h = {};
var file = new File(arguments[1]);
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (t[3] != gene) continue;
	if (h[t[0]] == null) h[t[0]] = [];
	h[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
}
file.close();

var s = {}, re = /(\d+)([MIDSHN])/g;
file = new File(arguments[2]);
while (file.readline(buf) >= 0) {
	var line = buf.toString();
	var m, t = line.split("\t");
	var x = h[t[2]];
	if (x == null) continue;

	var start = parseInt(t[3]) - 1, end = start;
	while ((m = re.exec(t[5])) != null) // parse CIGAR to get the end position
		if (m[2] == 'M' || m[2] == 'D')
			end += parseInt(m[1]);

	var max_ovlp = 0;
	for (var i = 0; i < x.length; ++i) {
		var max_left = x[i][0] > start? x[i][0] : start;
		var min_rght = x[i][1] < end  ? x[i][1] : end;
		max_ovlp = max_ovlp > min_rght - max_left? max_ovlp : min_rght - max_left;
	}

	var AS = null, XS = null;
	if ((m = /AS:i:(\d+)/.exec(line)) != null) AS = parseInt(m[1]);
	if ((m = /XS:i:(\d+)/.exec(line)) != null) XS = parseInt(m[1]);

	if (s[t[0]] == null) s[t[0]] = [];
	s[t[0]].push([AS, XS, max_ovlp]);
}
file.close();

buf.destroy();

for (var x in s) {
	var is_rejected = false, y = s[x];
	y.sort(function(a,b) {return b[0]-a[0]});
	for (var i = 0; i < y.length && y[i][0] == y[0][0]; ++i)
		if (y[0][2] < min_ovlp || y[i][0] == y[i][1])
			is_rejected = true;
	if (is_rejected) continue;
	print(x);
}
