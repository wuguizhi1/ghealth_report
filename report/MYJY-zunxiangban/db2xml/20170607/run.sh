rm -rf XML

perl listGenerateByProduct.pl -cus 常道康 -ser 简套系  -set 尊享版  -lan CN -out .

perl all2xml.pl -host 10.0.0.8 -port 27021 -db1 susceptibility -cl1 products -db2 susceptibility -cl2 prodata -lan CN -cus 常道康 -ser 简套系 -set 尊享版 -out XML

