make orthoPoly
make test
for i in `seq 1 10`
do
  ./tmp #100 | grep "Time:"
done
