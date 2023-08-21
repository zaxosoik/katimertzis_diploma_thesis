function Tinertia = Tinertia(Mrec,Ds,Ds_tonos,wmega) % ropi logo adraneias

Tinertia = Mrec*(-Ds_tonos).*(-Ds)*wmega.^2;

end

% ypologizetai to ena kommati tis Tinertia
% to kommati tou sintelesti me thita'' mpainei kateutheian sto main
% (Constant h Variable Inertia, pinakas Q)
% epilego an thelo ypologismo me Contant Inetria, Variable Inertia sto main