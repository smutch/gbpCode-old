#!/bin/csh 

set DIR_IN_ROOT=$1
set DIR_OUT_ROOT=$2
set NFILES=$3
set SNAP_NUM_LO=$4
set SNAP_NUM_HI=$5

set SNAP_NUM=$SNAP_NUM_LO
while ( $SNAP_NUM <= $SNAP_NUM_HI )
   if ( $SNAP_NUM > 99 ) then
      set SNAP_TXT = $SNAP_NUM
   else if ( $SNAP_NUM > 9 ) then
      set SNAP_TXT = "0"$SNAP_NUM
   else
      set SNAP_TXT = "00"$SNAP_NUM
   endif

   set I_SUFFIX=0
   while( $I_SUFFIX < 4)
      if ( $I_SUFFIX == 0) then
         set SUFFIX="catalog_groups_properties"
      else if ( $I_SUFFIX == 1) then
         set SUFFIX="catalog_groups_profiles"
      else if ( $I_SUFFIX == 2) then
         set SUFFIX="catalog_subgroups_properties"
      else if ( $I_SUFFIX == 3) then
         set SUFFIX="catalog_subgroups_profiles"
      endif
      set DIR_IN  = $DIR_IN_ROOT"_"$SNAP_TXT"."$SUFFIX
      set DIR_OUT = $DIR_OUT_ROOT"_"$SNAP_TXT"."$SUFFIX
      mv $DIR_IN $DIR_OUT
      set I_FILE=0
      while ( $I_FILE < $NFILES)
         set FILE_IN=$DIR_OUT"/"$DIR_IN_ROOT"_"$SNAP_TXT"."$SUFFIX"."$I_FILE
         set FILE_OUT=$DIR_OUT"/"$DIR_OUT_ROOT"_"$SNAP_TXT"."$SUFFIX"."$I_FILE
         mv $FILE_IN $FILE_OUT
         @ I_FILE = $I_FILE + 1
      end
      @ I_SUFFIX = $I_SUFFIX + 1
   end

   echo "Snapshot "$SNAP_TXT" Done."
   @ SNAP_NUM = $SNAP_NUM + 1
end

