# Run this script from easal-dev/processedLogs
# Takes two arguments:
#	$1: ouput csvfile
#	$2: Logs path
echo "ID,Min,1st_Qu,Median,Mean,3rd_Qu,Max,Std_Dev" > "${1}_part3.csv"
echo "ID,Min,1st_Qu,Median,Mean,3rd_Qu,Max,Std_Dev" > "${1}_part4.csv"
for i in ${2}/*; do 
  echo $i
  outname=`basename $i`
  sh ../scripts/extractLogs.sh -o=$outname $i/easal* 
  Rscript ../scripts/analyzeWaitTimesWrapper.R $outname $1
done 

