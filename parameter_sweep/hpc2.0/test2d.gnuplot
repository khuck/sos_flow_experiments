#set dgrid3d 30,30
set dgrid3d
set terminal postscript eps enhanced color font 'Helvetica,20'

set key outside horizontal top
set xlabel 'Doubles packed per publish'
set ylabel 'Latency (seconds)'
#set yrange [-0.25:0.25]
#set logscale y 2
# set the iso orientation
# set view 60,45
set xtics (1, 64, 128, 192, 256, 320, 384, 448, 512)
set boxwidth 2.0
set style fill solid 1.0 border 0
set style histogram errorbars gap 2 lw 1
set style data histograms
#set xtics rotate by -45
#set bars 0.5

set title '100 publish events per second'
set output '100-doubles-per-second.eps'
set yrange [-0.5:4.0]
plot 'latencies-listener.txt' u ($4-14.0):(($5==10000 && $3==1)?$9:1/0):(sqrt($11)) with yerrorbars title '1 ranks', \
     'latencies-listener.txt' u ($4-12.0):(($5==10000 && $3==2)?$9:1/0):(sqrt($11)) with yerrorbars title '2 ranks', \
     'latencies-listener.txt' u ($4-10.0):(($5==10000 && $3==4)?$9:1/0):(sqrt($11)) with yerrorbars title '4 ranks', \
     'latencies-listener.txt' u ($4-8.0):(($5==10000 && $3==6)?$9:1/0):(sqrt($11)) with yerrorbars title '6 ranks', \
     'latencies-listener.txt' u ($4-6.0):(($5==10000 && $3==8)?$9:1/0):(sqrt($11)) with yerrorbars title '8 ranks', \
     'latencies-listener.txt' u ($4-4.0):(($5==10000 && $3==10)?$9:1/0):(sqrt($11)) with yerrorbars title '10 ranks', \
     'latencies-listener.txt' u ($4-2.0):(($5==10000 && $3==12)?$9:1/0):(sqrt($11)) with yerrorbars title '12 ranks', \
     'latencies-listener.txt' u ($4):(($5==10000 && $3==14)?$9:1/0):(sqrt($11)) with yerrorbars title '14 ranks', \
     'latencies-listener.txt' u ($4+2.0):(($5==10000 && $3==16)?$9:1/0):(sqrt($11)) with yerrorbars title '16 ranks', \
     'latencies-listener.txt' u ($4+4.0):(($5==10000 && $3==18)?$9:1/0):(sqrt($11)) with yerrorbars title '18 ranks', \
     'latencies-listener.txt' u ($4+6.0):(($5==10000 && $3==20)?$9:1/0):(sqrt($11)) with yerrorbars title '20 ranks', \
     'latencies-listener.txt' u ($4+8.0):(($5==10000 && $3==22)?$9:1/0):(sqrt($11)) with yerrorbars title '22 ranks', \
     'latencies-listener.txt' u ($4+10.0):(($5==10000 && $3==24)?$9:1/0):(sqrt($11)) with yerrorbars title '24 ranks', \
     'latencies-listener.txt' u ($4+12.0):(($5==10000 && $3==26)?$9:1/0):(sqrt($11)) with yerrorbars title '26 ranks', \
     'latencies-listener.txt' u ($4+14.0):(($5==10000 && $3==28)?$9:1/0):(sqrt($11)) with yerrorbars title '28 ranks'

set title '10 publish events per second'
set output '10-doubles-per-second.eps'
set yrange [-0.10:0.25]
plot 'latencies-listener.txt' u ($4-14.0):(($5==100000 && $3==1)?$9:1/0):(sqrt($11)) with yerrorbars title '1 ranks', \
     'latencies-listener.txt' u ($4-12.0):(($5==100000 && $3==2)?$9:1/0):(sqrt($11)) with yerrorbars title '2 ranks', \
     'latencies-listener.txt' u ($4-10.0):(($5==100000 && $3==4)?$9:1/0):(sqrt($11)) with yerrorbars title '4 ranks', \
     'latencies-listener.txt' u ($4-8.0):(($5==100000 && $3==6)?$9:1/0):(sqrt($11)) with yerrorbars title '6 ranks', \
     'latencies-listener.txt' u ($4-6.0):(($5==100000 && $3==8)?$9:1/0):(sqrt($11)) with yerrorbars title '8 ranks', \
     'latencies-listener.txt' u ($4-4.0):(($5==100000 && $3==10)?$9:1/0):(sqrt($11)) with yerrorbars title '10 ranks', \
     'latencies-listener.txt' u ($4-2.0):(($5==100000 && $3==12)?$9:1/0):(sqrt($11)) with yerrorbars title '12 ranks', \
     'latencies-listener.txt' u ($4):(($5==100000 && $3==14)?$9:1/0):(sqrt($11)) with yerrorbars title '14 ranks', \
     'latencies-listener.txt' u ($4+2.0):(($5==100000 && $3==16)?$9:1/0):(sqrt($11)) with yerrorbars title '16 ranks', \
     'latencies-listener.txt' u ($4+4.0):(($5==100000 && $3==18)?$9:1/0):(sqrt($11)) with yerrorbars title '18 ranks', \
     'latencies-listener.txt' u ($4+6.0):(($5==100000 && $3==20)?$9:1/0):(sqrt($11)) with yerrorbars title '20 ranks', \
     'latencies-listener.txt' u ($4+8.0):(($5==100000 && $3==22)?$9:1/0):(sqrt($11)) with yerrorbars title '22 ranks', \
     'latencies-listener.txt' u ($4+10.0):(($5==100000 && $3==24)?$9:1/0):(sqrt($11)) with yerrorbars title '24 ranks', \
     'latencies-listener.txt' u ($4+12.0):(($5==100000 && $3==26)?$9:1/0):(sqrt($11)) with yerrorbars title '26 ranks', \
     'latencies-listener.txt' u ($4+14.0):(($5==100000 && $3==28)?$9:1/0):(sqrt($11)) with yerrorbars title '28 ranks'

set title '1 publish events per second'
set output '1-doubles-per-second.eps'
set yrange [-0.10:0.25]
plot 'latencies-listener.txt' u ($4-14.0):(($5==1000000 && $3==1)?$9:1/0):(sqrt($11)) with yerrorbars title '1 ranks', \
     'latencies-listener.txt' u ($4-12.0):(($5==1000000 && $3==2)?$9:1/0):(sqrt($11)) with yerrorbars title '2 ranks', \
     'latencies-listener.txt' u ($4-10.0):(($5==1000000 && $3==4)?$9:1/0):(sqrt($11)) with yerrorbars title '4 ranks', \
     'latencies-listener.txt' u ($4-8.0):(($5==1000000 && $3==6)?$9:1/0):(sqrt($11)) with yerrorbars title '6 ranks', \
     'latencies-listener.txt' u ($4-6.0):(($5==1000000 && $3==8)?$9:1/0):(sqrt($11)) with yerrorbars title '8 ranks', \
     'latencies-listener.txt' u ($4-4.0):(($5==1000000 && $3==10)?$9:1/0):(sqrt($11)) with yerrorbars title '10 ranks', \
     'latencies-listener.txt' u ($4-2.0):(($5==1000000 && $3==12)?$9:1/0):(sqrt($11)) with yerrorbars title '12 ranks', \
     'latencies-listener.txt' u ($4):(($5==1000000 && $3==14)?$9:1/0):(sqrt($11)) with yerrorbars title '14 ranks', \
     'latencies-listener.txt' u ($4+2.0):(($5==1000000 && $3==16)?$9:1/0):(sqrt($11)) with yerrorbars title '16 ranks', \
     'latencies-listener.txt' u ($4+4.0):(($5==1000000 && $3==18)?$9:1/0):(sqrt($11)) with yerrorbars title '18 ranks', \
     'latencies-listener.txt' u ($4+6.0):(($5==1000000 && $3==20)?$9:1/0):(sqrt($11)) with yerrorbars title '20 ranks', \
     'latencies-listener.txt' u ($4+8.0):(($5==1000000 && $3==22)?$9:1/0):(sqrt($11)) with yerrorbars title '22 ranks', \
     'latencies-listener.txt' u ($4+10.0):(($5==1000000 && $3==24)?$9:1/0):(sqrt($11)) with yerrorbars title '24 ranks', \
     'latencies-listener.txt' u ($4+12.0):(($5==1000000 && $3==26)?$9:1/0):(sqrt($11)) with yerrorbars title '26 ranks', \
     'latencies-listener.txt' u ($4+14.0):(($5==1000000 && $3==28)?$9:1/0):(sqrt($11)) with yerrorbars title '28 ranks'

set title '100 publish events per second'
set output '100-doubles-per-second-aggregator.eps'
set yrange [-0.5:4.0]
plot 'latencies-aggregator.txt' u ($4-14.0):(($5==10000 && $3==1)?$9:1/0):(sqrt($11)) with yerrorbars title '1 ranks', \
     'latencies-aggregator.txt' u ($4-12.0):(($5==10000 && $3==2)?$9:1/0):(sqrt($11)) with yerrorbars title '2 ranks', \
     'latencies-aggregator.txt' u ($4-10.0):(($5==10000 && $3==4)?$9:1/0):(sqrt($11)) with yerrorbars title '4 ranks', \
     'latencies-aggregator.txt' u ($4-8.0):(($5==10000 && $3==6)?$9:1/0):(sqrt($11)) with yerrorbars title '6 ranks', \
     'latencies-aggregator.txt' u ($4-6.0):(($5==10000 && $3==8)?$9:1/0):(sqrt($11)) with yerrorbars title '8 ranks', \
     'latencies-aggregator.txt' u ($4-4.0):(($5==10000 && $3==10)?$9:1/0):(sqrt($11)) with yerrorbars title '10 ranks', \
     'latencies-aggregator.txt' u ($4-2.0):(($5==10000 && $3==12)?$9:1/0):(sqrt($11)) with yerrorbars title '12 ranks', \
     'latencies-aggregator.txt' u ($4):(($5==10000 && $3==14)?$9:1/0):(sqrt($11)) with yerrorbars title '14 ranks', \
     'latencies-aggregator.txt' u ($4+2.0):(($5==10000 && $3==16)?$9:1/0):(sqrt($11)) with yerrorbars title '16 ranks', \
     'latencies-aggregator.txt' u ($4+4.0):(($5==10000 && $3==18)?$9:1/0):(sqrt($11)) with yerrorbars title '18 ranks', \
     'latencies-aggregator.txt' u ($4+6.0):(($5==10000 && $3==20)?$9:1/0):(sqrt($11)) with yerrorbars title '20 ranks', \
     'latencies-aggregator.txt' u ($4+8.0):(($5==10000 && $3==22)?$9:1/0):(sqrt($11)) with yerrorbars title '22 ranks', \
     'latencies-aggregator.txt' u ($4+10.0):(($5==10000 && $3==24)?$9:1/0):(sqrt($11)) with yerrorbars title '24 ranks', \
     'latencies-aggregator.txt' u ($4+12.0):(($5==10000 && $3==26)?$9:1/0):(sqrt($11)) with yerrorbars title '26 ranks', \
     'latencies-aggregator.txt' u ($4+14.0):(($5==10000 && $3==28)?$9:1/0):(sqrt($11)) with yerrorbars title '28 ranks'

set title '10 publish events per second'
set output '10-doubles-per-second-aggregator.eps'
set yrange [-0.10:0.25]
plot 'latencies-aggregator.txt' u ($4-14.0):(($5==100000 && $3==1)?$9:1/0):(sqrt($11)) with yerrorbars title '1 ranks', \
     'latencies-aggregator.txt' u ($4-12.0):(($5==100000 && $3==2)?$9:1/0):(sqrt($11)) with yerrorbars title '2 ranks', \
     'latencies-aggregator.txt' u ($4-10.0):(($5==100000 && $3==4)?$9:1/0):(sqrt($11)) with yerrorbars title '4 ranks', \
     'latencies-aggregator.txt' u ($4-8.0):(($5==100000 && $3==6)?$9:1/0):(sqrt($11)) with yerrorbars title '6 ranks', \
     'latencies-aggregator.txt' u ($4-6.0):(($5==100000 && $3==8)?$9:1/0):(sqrt($11)) with yerrorbars title '8 ranks', \
     'latencies-aggregator.txt' u ($4-4.0):(($5==100000 && $3==10)?$9:1/0):(sqrt($11)) with yerrorbars title '10 ranks', \
     'latencies-aggregator.txt' u ($4-2.0):(($5==100000 && $3==12)?$9:1/0):(sqrt($11)) with yerrorbars title '12 ranks', \
     'latencies-aggregator.txt' u ($4):(($5==100000 && $3==14)?$9:1/0):(sqrt($11)) with yerrorbars title '14 ranks', \
     'latencies-aggregator.txt' u ($4+2.0):(($5==100000 && $3==16)?$9:1/0):(sqrt($11)) with yerrorbars title '16 ranks', \
     'latencies-aggregator.txt' u ($4+4.0):(($5==100000 && $3==18)?$9:1/0):(sqrt($11)) with yerrorbars title '18 ranks', \
     'latencies-aggregator.txt' u ($4+6.0):(($5==100000 && $3==20)?$9:1/0):(sqrt($11)) with yerrorbars title '20 ranks', \
     'latencies-aggregator.txt' u ($4+8.0):(($5==100000 && $3==22)?$9:1/0):(sqrt($11)) with yerrorbars title '22 ranks', \
     'latencies-aggregator.txt' u ($4+10.0):(($5==100000 && $3==24)?$9:1/0):(sqrt($11)) with yerrorbars title '24 ranks', \
     'latencies-aggregator.txt' u ($4+12.0):(($5==100000 && $3==26)?$9:1/0):(sqrt($11)) with yerrorbars title '26 ranks', \
     'latencies-aggregator.txt' u ($4+14.0):(($5==100000 && $3==28)?$9:1/0):(sqrt($11)) with yerrorbars title '28 ranks'

set title '1 publish events per second'
set output '1-doubles-per-second-aggregator.eps'
set yrange [-0.10:0.25]
plot 'latencies-aggregator.txt' u ($4-14.0):(($5==1000000 && $3==1)?$9:1/0):(sqrt($11)) with yerrorbars title '1 ranks', \
     'latencies-aggregator.txt' u ($4-12.0):(($5==1000000 && $3==2)?$9:1/0):(sqrt($11)) with yerrorbars title '2 ranks', \
     'latencies-aggregator.txt' u ($4-10.0):(($5==1000000 && $3==4)?$9:1/0):(sqrt($11)) with yerrorbars title '4 ranks', \
     'latencies-aggregator.txt' u ($4-8.0):(($5==1000000 && $3==6)?$9:1/0):(sqrt($11)) with yerrorbars title '6 ranks', \
     'latencies-aggregator.txt' u ($4-6.0):(($5==1000000 && $3==8)?$9:1/0):(sqrt($11)) with yerrorbars title '8 ranks', \
     'latencies-aggregator.txt' u ($4-4.0):(($5==1000000 && $3==10)?$9:1/0):(sqrt($11)) with yerrorbars title '10 ranks', \
     'latencies-aggregator.txt' u ($4-2.0):(($5==1000000 && $3==12)?$9:1/0):(sqrt($11)) with yerrorbars title '12 ranks', \
     'latencies-aggregator.txt' u ($4):(($5==1000000 && $3==14)?$9:1/0):(sqrt($11)) with yerrorbars title '14 ranks', \
     'latencies-aggregator.txt' u ($4+2.0):(($5==1000000 && $3==16)?$9:1/0):(sqrt($11)) with yerrorbars title '16 ranks', \
     'latencies-aggregator.txt' u ($4+4.0):(($5==1000000 && $3==18)?$9:1/0):(sqrt($11)) with yerrorbars title '18 ranks', \
     'latencies-aggregator.txt' u ($4+6.0):(($5==1000000 && $3==20)?$9:1/0):(sqrt($11)) with yerrorbars title '20 ranks', \
     'latencies-aggregator.txt' u ($4+8.0):(($5==1000000 && $3==22)?$9:1/0):(sqrt($11)) with yerrorbars title '22 ranks', \
     'latencies-aggregator.txt' u ($4+10.0):(($5==1000000 && $3==24)?$9:1/0):(sqrt($11)) with yerrorbars title '24 ranks', \
     'latencies-aggregator.txt' u ($4+12.0):(($5==1000000 && $3==26)?$9:1/0):(sqrt($11)) with yerrorbars title '26 ranks', \
     'latencies-aggregator.txt' u ($4+14.0):(($5==1000000 && $3==28)?$9:1/0):(sqrt($11)) with yerrorbars title '28 ranks'

