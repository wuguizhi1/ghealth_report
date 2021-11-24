perl get_result_xml.pl -data tptest/ -barcode tptest001 -out tptest/tptest001_raw && perl meta_make_pdf_v2.pl -xml tptest/tptest001_raw/tptest001.xml -barcode tptest001 -od tptest/tptest001_raw/pdf
perl get_result_xml.pl -data tptest/ -barcode tptest001 -l l2 -out tptest/tptest001_l2 && perl meta_make_pdf_v2.pl -xml tptest/tptest001_l2/tptest001.xml -barcode tptest001 -od tptest/tptest001_l2/pdf
perl get_result_xml.pl -data tptest/ -barcode tptest001 -l l3 -out tptest/tptest001_l3 && perl meta_make_pdf_v2.pl -xml tptest/tptest001_l3/tptest001.xml -barcode tptest001 -od tptest/tptest001_l3/pdf
