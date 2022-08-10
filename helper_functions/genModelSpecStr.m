function [modelSpecStr] = genModelSpecStr(agentType,actionSelectionMethod,...
                                          utilityFunc1,utilityFunc2,...
                                          initializationMethod,forgettingType,...
                                          fullANS,noAbortANS)
modelSpecStr = [agentType '_' actionSelectionMethod];
if (~strcmp(utilityFunc1,''))
    modelSpecStr = [modelSpecStr '_' utilityFunc1];
end
if (~strcmp(utilityFunc2,''))
    modelSpecStr = [modelSpecStr '_' utilityFunc2];
end
modelSpecStr = [modelSpecStr '_initializationMethod_' initializationMethod '_forgettingType_' forgettingType];
if (fullANS)
    modelSpecStr = [modelSpecStr '_fullANS'];
end
if (noAbortANS)
    modelSpecStr = [modelSpecStr '_noAbortANS'];
end

end

