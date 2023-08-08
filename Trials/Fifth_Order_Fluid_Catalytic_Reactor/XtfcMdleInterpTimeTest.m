function ric = BsTimeTest(timeStore, xiStore, SfStore, inputStruct, mdle3DMats)

CE = TFC_MDRE_Interp(0.4711, timeStore, xiStore, SfStore, inputStruct);
CE3 = reshape(CE', 5, 5, 1);
ric = (mdle3DMats.Sm(:,:,1) + inv(CE3(:,:)));

end