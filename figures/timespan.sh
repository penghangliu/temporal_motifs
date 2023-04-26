#- take each output file and create as many distributions as the number of motif variants
#- to handle the overlapping events (all dc30s will appear in dc60), you can put all (dc values) in the same chart and start plotting with the largest dc -- smaller dc values will be on top for overlapping events

motif=$1
data=$2

grep span $data|grep $motif|awk '{print $2}'|sort -n >hist

font=20
echo "set terminal postscript landscape enhanced color font 'Helvetica,$font'" >> tmp
echo "set title \"$motif\"" >>tmp

name=`echo $data | cut -c1-6`
echo "set output '"$name"-"$motif"_spanHist.ps'" >> tmp

echo "set ylabel 'counts'" >>tmp

echo "set xlabel 'end-to-end timespan'" >> tmp

echo "set boxwidth 1" >>tmp
echo "set style fill solid 1.0 border rgb 'blue'" >>tmp
echo "plot 'hist' using 1 bins=100 with boxes lc rgb 'blue' notitle" >>tmp

gnuplot tmp
rm tmp
ps2pdf $name"-"$motif"_spanHist.ps"
rm $name"-"$motif"_spanHist.ps"

rm hist
