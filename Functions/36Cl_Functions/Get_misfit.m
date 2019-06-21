function RMSw = Get_misfit(x0,data_mc)
% Function to get the RMSw misfit of a given model

% get theoretical 36Cl for each sample                              
N_36 = Model_direct_36Facet(x0,data_mc);
% get the RMSw                              
RMSw = sum((((data_mc.dataset(1,:)-N_36).^2).^.5)./data_mc.dataset(2,:));
end