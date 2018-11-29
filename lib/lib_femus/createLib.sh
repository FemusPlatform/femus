#
# SCRIPT FOR FEMUS LIB
#

make

mv libfemus.so $FEMUS_DIR/lib/
echo "Linking femus lib in $FEMUS_DIR/lib/ directory"
