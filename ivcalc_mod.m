function iv = ivcalc_mod(obj,Ind)
   [FreqTab, nanFreqTab]   = obj.BinContainers{Ind}.getFrequencyTable();

   [iv,ivrows] = InformationValueByFreq_mod(FreqTab, nanFreqTab);
end
    
    