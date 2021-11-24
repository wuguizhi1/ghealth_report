rm -rf XML

perl listGenerateByProduct.pl -cus GI-T客户 -ser 单项系  -set 体重管理  -lan CN -out .

perl all2xml.pl -host 10.0.0.204 -port 27021 -db1 susceptibility -cl1 products -db2 susceptibility -cl2 prodata -lan CN -cus GI-T客户 -ser 单项系 -set 体重管理 -out XML

