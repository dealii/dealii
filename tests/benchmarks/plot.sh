


echo "set terminal x11 persist"
echo "set log y"
echo "plot \\"
n=1
for i in `cat names.$1`
do
  n=`expr $n "+" 1`
  echo "'datatable.$1' using 1:$n title '$i' w lp,\\";
done

echo "0"

echo "pause -1"

