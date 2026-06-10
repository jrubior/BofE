function [St_old] = restriction_kilianmurphy(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info)
% checking sign restrictions


    %BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

    %structpara   = fo_inv(BSigmaQnew);

    %irfpara      = fo_str2irfs(structpara);


    hSigmadraw = chol(Sigmadraw);

    L0 = hSigmadraw'*Qdraw;%reshape(irfpara(1:info.nvar*info.nvar,:),info.nvar,info.nvar);


    supply_shock = (L0(1,1)<0)*(L0(2,1)<0)*(L0(3,1)>0);
    flow_demand_shock = (L0(1,2)>0)*(L0(2,2)>0)*(L0(3,2)>0);
    spec_demand_shock = (L0(1,3)>0)*(L0(2,3)<0)*(L0(3,3)>0)*(L0(4,3)>0);


    elasticity=L0(1,3)./L0(3,3);   %supply elasticity in response to speculative demand shock
    ADelas    =L0(1,2)./L0(3,2);   %supply elasticity in response to flow demand shock


    elast_restrictions = (ADelas<=.0258)*(elasticity<=.0258);%(elasticity<=.0258)*(ADelas<=.0258);



    % dynamic impulse responses

    A0         = hSigmadraw\Qdraw;
    Aplus      = Bdraw*A0;
    for l=1:info.nlag-1
        info.A{l} = Aplus((l-1)*info.nvar+info.nex+1:l*info.nvar+info.nex,1:end);
        info.bigF((l-1)*info.nvar+1:l*info.nvar,1:info.nvar)=info.A{l}/A0;
    end
    A{info.nlag} = Aplus((info.nlag-1)*info.nvar+1:info.nlag*info.nvar,1:end);
    info.bigF((info.nlag-1)*info.nvar+1:info.nlag*info.nvar,:)=[A{info.nlag}/A0 info.extraF];


    
    activity_index = zeros(12,1);
    oil_price      = zeros(12,1);

    for h=0:1:11

        tmp = (info.J'*((info.bigF')^h)*info.J)*L0;
        
        activity_index(h+1) =tmp(2,1);
        oil_price(h+1)      =tmp(3,1);

    end


    label= supply_shock*flow_demand_shock*spec_demand_shock*elast_restrictions*all(activity_index<0)*all(oil_price>0);


    


    St_old = label;








