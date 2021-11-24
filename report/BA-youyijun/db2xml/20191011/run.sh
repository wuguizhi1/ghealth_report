rm -rf XML

perl listGenerateByProduct.pl -cus 常道康 -ser 博奥 -set 有益菌 -lan CN -out .

perl all2xml.pl -host 10.0.0.204 -port 27021 -db1 susceptibility -cl1 products -db2 susceptibility -cl2 prodata -lan CN -cus 常道康 -ser 博奥 -set 有益菌 -out XML

