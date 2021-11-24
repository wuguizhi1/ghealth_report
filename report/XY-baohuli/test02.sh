perl get_result_xml.pl -data test/ -barcode tptest001 -out test/tptest001_raw && perl meta_make_pdf_v2.pl -xml test/tptest001_raw/tptest001.xml -barcode tptest001 -od test/tptest001_raw/pdf
perl get_result_xml.pl -data test/ -barcode tptest002 -out test/tptest002_raw && perl meta_make_pdf_v2.pl -xml test/tptest002_raw/tptest002.xml -barcode tptest002 -od test/tptest002_raw/pdf
perl get_result_xml.pl -data test/ -barcode tptest003 -out test/tptest003_raw && perl meta_make_pdf_v2.pl -xml test/tptest003_raw/tptest003.xml -barcode tptest003 -od test/tptest003_raw/pdf
