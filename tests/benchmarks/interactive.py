import textwrap
import os

begin = \
"""
<!DOCTYPE HTML>
<html>
	<head>
		<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
		<title>Highcharts Example</title>

		<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.8.2/jquery.min.js"></script>

<script src="http://code.highcharts.com/highcharts.js"></script>
<script src="http://code.highcharts.com/modules/exporting.js"></script>

<div id="container" style="min-width: 400px; height: 400px; margin: 0 auto"></div>
<script>
$(function () {
    var chart;
    $(document).ready(function() {
        chart = new Highcharts.Chart({
            chart: {
                renderTo: 'container',
                type: 'line',
                marginRight: 250,
                marginBottom: 25,
                zoomType: 'x'
            },
            title: {
                text: 'regression timings',
                x: -20 //center
            },
            xAxis: {
            },
            yAxis: {
                title: {
                    text: 'slowdown (%)'
                },
                plotLines: [{
                    value: 0,
                    width: 1,
                    color: '#808080'
                }]
            },
            tooltip: {
                formatter: function() {
                        return '<b>'+ this.series.name +'</b><br/>'+
                        this.y +'s' + ', rev ' + this.x;
                }
            },
            legend: {
                layout: 'vertical',
                align: 'right',
                verticalAlign: 'top',
                x: -10,
                y: 100,
                borderWidth: 0
            },

            series: [
"""

end = \
"""
]
        });
setTimeout(function(){
chart.xAxis[0].setExtremes(28000,chart.xAxis[0].getExtremes().dataMax); 
},1000);
    });
});    


</script>
	</head>
	<body>
deal.II performance benchmarks, see 
<a href="http://www.dealii.org/testsuite.html">http://www.dealii.org/testsuite.html</a><br>

<div id="container" style="min-width: 400px; height: 400px; margin: 0 auto"></div>

<img src="baseline.png"/>

<img src="test_assembly.png"/>

<img src="step-22.png"/>

<img src="test_poisson.png"/>

<img src="tablehandler.png"/>


	</body>
</html>
"""



list = os.listdir(".")

print begin

first = 1
for fname in list:
    if (fname.startswith("names.")):
	testname = fname[6:]
	names = open(fname).readlines()
	data = open("datatable."+testname).readlines()
	idx = 0
	for name in names:
	    if first == 1:
		first = 0
	    else:
		print ","
	    idx = idx+1
	    print "{ name: '%s - %s', data: [" % (testname,name[:-1])
	    i=0
#((\$$col-$baseline)/$baseline*100.0)
	    baseline = -1
	    for l in data:
	        if (len(l.strip())<1):
		    continue
		
		lnumbers = [float(x) for x in l.split()]
		if len(lnumbers)<=1:
		    continue;

		if baseline<0:
		    baseline=lnumbers[idx]

		if lnumbers[0]==0:
		    continue;
		
		if (i==1):
		    print(","),
		i=1;
		    
		print "[%d,%f]" % (lnumbers[0], (lnumbers[idx]-baseline)/baseline*100.0)
	    print "]}\n"
       



print end