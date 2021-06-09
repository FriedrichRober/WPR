#
# WPR: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "WPR" );

TestDirectory(DirectoriesPackageLibrary( "WPR", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
