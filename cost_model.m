function [costSys] = cost_model( W_empty,Velocity_max,Thrust_SLS,Temp_Turbine,quantity,year,month )
% cost_model is a function created to calculate cost parameter for an
% aircraft system. 
%-------------------------------INPUT--------------------------------------%
% W_empty --> Empty weight of the system
% Velocity_max --> Maximum velocity for the system
% Thrust_SLS --> Max thrust for engine at sea level
% numOfEngine --> Total number of engines onboard the system
% quantity --> number of aircraft being produced
% year --> year of producation
% month --> month of producation, if month = 0 then it takes a running
% average of the inflation for the year
%inflation = cpi_search(year,month);
inflation = 1;
% hourly rates
HR_Tooling = 2.883*year-5666;
HR_Engineering = 2.576*year-5058;
HR_Manufacturing = 2.316*year-4552;
HR_QualityControl = 2.6*year-5112;
% tooling hours
H_Tooling = 5.99.*W_empty^0.777.*Velocity_max.^0.696.*quantity.^0.263;
% Airframe engineering hours
H_Engineering = 4.86.*W_empty.^0.777.*Velocity_max.^0.894.*quantity.^0.163;
% Manufacturing labor hours
H_Manufacturing = 7.37.*W_empty.^0.82.*Velocity_max.^0.484.*quantity.^0.641;
% Quality control hours for cargo and transport aircraft
 H_QualityControl = 0.076.*H_Manufacturing;
% H_QualityControl = 0.13.*H_Manufacturing;

costSys.COST_TOOLING = H_Tooling.*HR_Tooling;
costSys.COST_ENGINEERING = H_Engineering.*HR_Engineering;
costSys.COST_MANUFACTURING = H_Manufacturing.*HR_Manufacturing;
costSys.COST_QUALITYCONTROL = H_QualityControl.*HR_QualityControl;


% Development support cost 
costSys.COST_DEVELOPMENTSUPPORT = inflation.*66.*W_empty.^0.63.*Velocity_max.^1.3;
% Flight Operations cost
costSys.COST_FLIGHTTESTOPERATIONS = inflation.*687.*W_empty.^0.325.*Velocity_max.^0.822.*quantity.^1.21;
% Manufacturing material and equipment cost
costSys.COST_MANUFACTURINGMATANDEQUIP = inflation.*16.39.*W_empty.^0.921.*Velocity_max.^0.621.*quantity.^0.799;
% Producation engine unit cost in 1998 dollars
% Total cost of the system
costSys.COST = costSys.COST_ENGINEERING + costSys.COST_DEVELOPMENTSUPPORT +...
    costSys.COST_FLIGHTTESTOPERATIONS + costSys.COST_TOOLING +...
		costSys.COST_MANUFACTURING + costSys.COST_QUALITYCONTROL...
        + costSys.COST_MANUFACTURINGMATANDEQUIP;
%    costSys = [COST_TOOLING , COST_ENGINEERING ,COST_MANUFACTURING , COST_QUALITYCONTROL ,...
%        COST_DEVELOPMENTSUPPORT , COST_FLIGHTTESTOPERATIONS , COST_MANUFACTURINGMATANDEQUIP, COST];

    
  



end
