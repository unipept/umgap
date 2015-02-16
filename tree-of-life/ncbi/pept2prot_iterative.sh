i=1
for file in results/fastafiles/*
do
  echo "Processing file $i"
  filename=$(basename $file)
  cat $file | unipept pept2prot --extra > "results/pept2prot_fastafiles/$filename"
  i=$((i+1))
done
