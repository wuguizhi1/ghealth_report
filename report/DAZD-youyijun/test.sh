perl get_result_xml.pl -data test/ -barcode GI00001 -out test/GI00001_raw && perl meta_make_pdf_v2.pl -xml test/GI00001_raw/GI00001.xml -barcode GI00001 -od test/GI00001_raw/pdf
perl get_result_xml.pl -data test/ -barcode GI00001 -l l2 -out test/GI00001_l2 && perl meta_make_pdf_v2.pl -xml test/GI00001_l2/GI00001.xml -barcode GI00001 -od test/GI00001_l2/pdf
perl get_result_xml.pl -data test/ -barcode GI00001 -l l3 -out test/GI00001_l3 && perl meta_make_pdf_v2.pl -xml test/GI00001_l3/GI00001.xml -barcode GI00001 -od test/GI00001_l3/pdf
