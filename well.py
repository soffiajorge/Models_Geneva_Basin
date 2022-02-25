from paraview.simple import *
well_list = ['P01','P02','P03','P04','P05','P06','P07','P08','P09','P10','P11','P12','P13','P14','I01','I02','I03','I04','I05','I06','I07','I08','I09','I10','I10'];
well_type = ['PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','PRODUCER','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR','INJECTOR'];
well_x = [38,21,44,31,33,19,15,36,46,50,27,65,61,57,49,31,48,59,55,36,33,29,24,48,42];
well_y = [36,36,43,27,18,30,40,42,23,18,41,23,35,23,23,19,34,17,30,28,39,41,28,11,18];
DX = 100;
DY = 100;
X0 = 350858.5624
Y0 = 7513812.6952
well_cyl = Cylinder();
SetProperties(well_cyl,Height=1000,Radius=30);
for idx, val in enumerate(well_list):
	t = Transform(well_cyl);
	t.Transform.Translate=[X0 + well_x[idx]*DX, Y0 + well_y[idx]*DY,3100];
	t.Transform.Rotate = [90,0,0];
	dp = GetDisplayProperties(t);
	if (well_type[idx] == 'PRODUCER'):
		dp.DiffuseColor=[1,0,0];
	else:	
		dp.DiffuseColor=[0,0,1];
	Show(t);
Render();

	
	





