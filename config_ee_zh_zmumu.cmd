!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	Main:numberOfEvents = 10000         ! number of events to generate
	
	Beams:idA = 11                   ! first beam, e- = -11
	Beams:idB = -11                  ! second beam, e+ = 11
	Beams:eCM = 240.                 ! CM energy of collision
	
	! Higgsstrahlung process
	HiggsSM:ffbar2HZ = on
	
	! 5) Force the Z decays to muons
	23:onMode = off
	23:onIfAny = 13 -13
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
