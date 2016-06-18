function composelocid_test()

locid1 = [63450002;120450002];

[statecode,countycode,siteid] = decomposelocid(locid1);

locid2 = composelocid(statecode,countycode,siteid);

disp(['input1: ',num2str(locid1(1)),' output1: ',num2str(locid2(1)),...
      ' input2: ',num2str(locid1(2)),' output2: ',num2str(locid2(2))]);
