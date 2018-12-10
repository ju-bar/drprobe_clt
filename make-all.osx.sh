# LINUX version of making a new Dr. Probe - Command-Line Tool distribution
echo "Dr. Probe (CLT) - re-building ..."
# cleaning previous object code
echo "- cleaning previous object code"
make clean -C celslc
make clean -C msa
make clean -C wavimg
# archiving current source code
echo "- archiving current sources"
./compress-drprobe_clt-src.sh
# making new object code from current sources
echo "- making object code from current sources"
make -C celslc
make -C msa
make -C wavimg
# updating the binary folder
echo "- updateing the binary folder"
./update-bin.osx.sh
# archiving the binary folder
echo "- archiving the binary folder"
./compress-drprobe_clt-bin.osx.sh
echo "Done."
