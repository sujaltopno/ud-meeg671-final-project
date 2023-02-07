%Function for computing a values
function a = a_values(qf)

qi = deg2rad([58.2686 75.3224 11.7968 45.9029 -22.1081 -31.2831 -42.3712]);
a = zeros(7,2);

for i = 1:7
    a2 = (-3*(qi(i)-qf(i)))/25;
    a3 = (2*(qi(i)-qf(i)))/125;
    a(i,1) = a3;
    a(i,2) = a2;
end
           