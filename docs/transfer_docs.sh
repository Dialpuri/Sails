output_file=../../s.tree

echo """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<!DOCTYPE instance-profile
        SYSTEM \"https://resources.jetbrains.com/writerside/1.0/product-profile.dtd\">

<instance-profile id=\"s\"
                 name=\"Sails\"
                 start-page=\"Sails.md\">

  <toc-element topic=\"Sails.md\">

""" > $output_file

dir_path="./"
file_string=""
for file in $(ls -p $dir_path | grep -v /)
do
    if [ $file != "Sails.md" ]
    then
    # Append the file name to the string
      file_string+="\t<toc-element topic=\"$file\"/>\n"
    fi
done

echo $file_string >> $output_file

echo """
    </toc-element>
    <toc-element topic=\"Installation.topic\"/>
    <toc-element topic=\"Command-Line-Usage.topic\"/>
</instance-profile>
""" >> $output_file