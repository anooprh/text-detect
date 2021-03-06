function [ final ] = detecttext( imName )
%DETECTTEXT Summary of this function goes here
%   Detailed explanation goes here

image = imread(imName);
swtMap1 = swt(image,1);
[swtLabel1, numCC1] = swtlabel(swtMap1);
final1 = extractletters(swtMap1, swtLabel1, numCC1);
swtMap2 = swt(image,-1);
[swtLabel2, numCC2] = swtlabel(swtMap2);
final2 = extractletters(swtMap2, swtLabel2, numCC2);
final = final1 | final2;

end

