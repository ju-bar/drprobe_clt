echo "Updating the drprobe_clt source archive"
rm "drprobe_clt_src.tar.gz"
tar cvzf "drprobe_clt_src.tar.gz" "BuildCell" "CellMuncher_v2.0" "common" "celslc" "msa" "wavimg" "test" "README.md" "LICENSE"
echo Done.
