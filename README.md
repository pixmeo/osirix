# Fork of Osirix DICOM Viewer

This fork builds a 64 bit version of the viewer. To do this, Papyrus and DICOM print support have been removed, and the Kakadu JPEG 2000 library has been replaced with the OpenJPEG 2.0 library.

This build ls linked with ITK 4.4 and VTK 6.1.  

Currently, the 3D rendering is broken.

To build, clone this repo, then:

clone the osirixBinaries repo, rename to Binaries,
and place in the osirix folder.

Now you can build the "Osirix" target with XCode 5.

Enjoy!!
