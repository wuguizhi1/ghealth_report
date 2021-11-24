rm -rf XML

perl listGenerateByProduct.pl -cus GI-S客户 -ser 肠道微生态  -set 基础版  -lan CN -out .

perl all2xml.pl -host 10.0.0.8 -port 27021 -db1 susceptibility -cl1 products -db2 susceptibility -cl2 prodata -lan CN -cus GI-S客户 -ser 肠道微生态 -set 基础版 -out XML

