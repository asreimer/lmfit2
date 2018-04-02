The test.rawacf file was made by taking the first 144844 bytes of the "20120614.0401.00.sas.rawacf" file using the command:

head -c 144844 20120614.0401.00.sas.rawacf > test.rawacf

This corresponds to 4 rawacf records.

Future testing code should run make_lmfit2 using the test.rawacf file and then compare the output to the expected.lmfit2 file. It is expected that the origin.time and origin.command fields in the new lmfit2 and expected.lmfit2 files should differ.