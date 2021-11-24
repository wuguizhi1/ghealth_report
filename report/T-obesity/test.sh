rm -r test/Demo_*
perl get_result_xml.pl -data test/Demo/ -barcode Demo -out test/Demo_raw && perl meta_make_pdf_v2.pl -xml test/Demo_raw/Demo.xml -od test/Demo_raw -barcode Demo
perl get_result_xml.pl -data test/Demo/ -barcode Demo -out test/Demo_l1 -l l1 && perl meta_make_pdf_v2.pl -xml test/Demo_l1/Demo.xml -od test/Demo_l1 -barcode Demo
perl get_result_xml.pl -data test/Demo/ -barcode Demo -out test/Demo_l2 -l l2 && perl meta_make_pdf_v2.pl -xml test/Demo_l2/Demo.xml -od test/Demo_l2 -barcode Demo
perl get_result_xml.pl -data test/Demo/ -barcode Demo -out test/Demo_l3 -l l3 && perl meta_make_pdf_v2.pl -xml test/Demo_l3/Demo.xml -od test/Demo_l3 -barcode Demo
perl get_result_xml.pl -data test/Demo/ -barcode Demo -out test/Demo_l4 -l l4 && perl meta_make_pdf_v2.pl -xml test/Demo_l4/Demo.xml -od test/Demo_l4 -barcode Demo
perl get_result_xml.pl -data test/Demo/ -barcode Demo -out test/Demo_l5 -l l5 && perl meta_make_pdf_v2.pl -xml test/Demo_l5/Demo.xml -od test/Demo_l5 -barcode Demo
