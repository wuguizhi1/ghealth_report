perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./0005550074/0005550074.result.json -prefix 0005550074 -o ./0005550074/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./0005550294/0005550294.result.json -prefix 0005550294 -o ./0005550294/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./0005550331/0005550331.result.json -prefix 0005550331 -o ./0005550331/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./0005550499/0005550499.result.json -prefix 0005550499 -o ./0005550499/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./0005550514/0005550514.result.json -prefix 0005550514 -o ./0005550514/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./GI01008812/GI01008812.result.json -prefix GI01008812 -o ./GI01008812/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./GI01009152/GI01009152.result.json -prefix GI01009152 -o ./GI01009152/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./GI01009703/GI01009703.result.json -prefix GI01009703 -o ./GI01009703/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./GI01009710/GI01009710.result.json -prefix GI01009710 -o ./GI01009710/ 
perl /data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip/bin/../scripts/extract_sample.pl -json ./KY03006749/KY03006749.result.json -prefix KY03006749 -o ./KY03006749/ 

perl ../calculate_meta_level_per.pl -infile ./0005550074/0005550074.extract.xls -outfile ./0005550074/0005550074.level.xls
perl ../calculate_meta_level_per.pl -infile ./0005550294/0005550294.extract.xls -outfile ./0005550294/0005550294.level.xls
perl ../calculate_meta_level_per.pl -infile ./0005550331/0005550331.extract.xls -outfile ./0005550331/0005550331.level.xls
perl ../calculate_meta_level_per.pl -infile ./0005550499/0005550499.extract.xls -outfile ./0005550499/0005550499.level.xls
perl ../calculate_meta_level_per.pl -infile ./0005550514/0005550514.extract.xls -outfile ./0005550514/0005550514.level.xls
perl ../calculate_meta_level_per.pl -infile ./GI01008812/GI01008812.extract.xls -outfile ./GI01008812/GI01008812.level.xls
perl ../calculate_meta_level_per.pl -infile ./GI01009152/GI01009152.extract.xls -outfile ./GI01009152/GI01009152.level.xls
perl ../calculate_meta_level_per.pl -infile ./GI01009703/GI01009703.extract.xls -outfile ./GI01009703/GI01009703.level.xls
perl ../calculate_meta_level_per.pl -infile ./GI01009710/GI01009710.extract.xls -outfile ./GI01009710/GI01009710.level.xls
perl ../calculate_meta_level_per.pl -infile ./KY03006749/KY03006749.extract.xls -outfile ./KY03006749/KY03006749.level.xls

perl ../level_xml.v1.pl -data ./0005550074/ -level_per ./0005550074/0005550074.level.xls -out ./0005550074 -barcode 0005550074
perl ../level_xml.v1.pl -data ./0005550294/ -level_per ./0005550294/0005550294.level.xls -out ./0005550294 -barcode 0005550294
perl ../level_xml.v1.pl -data ./0005550331/ -level_per ./0005550331/0005550331.level.xls -out ./0005550331 -barcode 0005550331
perl ../level_xml.v1.pl -data ./0005550499/ -level_per ./0005550499/0005550499.level.xls -out ./0005550499 -barcode 0005550499
perl ../level_xml.v1.pl -data ./0005550514/ -level_per ./0005550514/0005550514.level.xls -out ./0005550514 -barcode 0005550514
perl ../level_xml.v1.pl -data ./GI01008812/ -level_per ./GI01008812/GI01008812.level.xls -out ./GI01008812 -barcode GI01008812
perl ../level_xml.v1.pl -data ./GI01009152/ -level_per ./GI01009152/GI01009152.level.xls -out ./GI01009152 -barcode GI01009152
perl ../level_xml.v1.pl -data ./GI01009703/ -level_per ./GI01009703/GI01009703.level.xls -out ./GI01009703 -barcode GI01009703
perl ../level_xml.v1.pl -data ./GI01009710/ -level_per ./GI01009710/GI01009710.level.xls -out ./GI01009710 -barcode GI01009710
perl ../level_xml.v1.pl -data ./KY03006749/ -level_per ./KY03006749/KY03006749.level.xls -out ./KY03006749 -barcode KY03006749
