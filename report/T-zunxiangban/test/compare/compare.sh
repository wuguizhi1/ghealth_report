rm -rf compare_report
perl ../../get_result_xml.pl -data KYCD2523 -barcode KYCD2523 -historydir ./ -out compare_report
perl ../../meta_make_pdf_v2.pl -xml compare_report/KYCD2523.xml -barcode KYCD2523 -od compare_report
