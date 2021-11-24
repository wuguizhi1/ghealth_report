rm -rf XML

perl listGenerateByProduct.pl -cus 常道康 -ser 单项系  -set 肠道保护力  -lan CN -out .

perl all2xml.pl -host 10.0.0.8 -port 27021 -db1 susceptibility -cl1 products -db2 susceptibility -cl2 prodata -lan CN -cus 常道康 -ser 单项系 -set 肠道保护力 -out XML

