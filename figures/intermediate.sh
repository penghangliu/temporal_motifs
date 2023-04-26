
motif=$1
data=$2

grep perc $data|grep $motif|awk '{print $3}'|sort -n|uniq -c >hist2
grep perc $data|grep $motif|awk '{print $4}'|sort -n|uniq -c >hist3

font=32
echo "set terminal postscript landscape enhanced color font 'Helvetica,$font'" >> tmp
echo "set title \"$motif\"" >>tmp

name=`echo $data | cut -c1-6`
echo "set output '"$name"-"$motif"_percHist.ps'" >> tmp
#echo "set output '"$motif"_percHist.ps'" >> tmp

echo "set ylabel 'counts'" >>tmp

#echo "set xlabel 'second event occurrence percentile'" >> tmp
echo "set xlabel 'second/third event occurrence perc.'" >> tmp
echo "set xrange [0:100]" >> tmp
echo "set xtics 10 font 'Helvetica, $font'" >> tmp

echo "set boxwidth 1" >>tmp
echo "set style fill solid 1.0 border" >>tmp

echo "plot 'hist2' using 2:1 with boxes fill solid 1.0 lc rgb 'blue' notitle,\
    'hist3' using 2:1 with boxes fill empty lc rgb 'red' notitle" >>tmp


gnuplot tmp
rm tmp
ps2pdf $name"-"$motif"_percHist.ps"
rm $name"-"$motif"_percHist.ps"

