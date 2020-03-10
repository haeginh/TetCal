# TetCal (VoxelVersion 이랑 코드 최대한 맞추는 것을 목표로 함)
Geant4 code to calculate radiation doses with tetrahedral mesh phanotms

* 아직 구현 안된 부분

surface source

* 구현된 부분

External, Internal, Hadronic, dose 정리, effective dose 정리

* dose 정리 관련

[팬텀이름].dose 라는 파일 내에 아래 형식에 맞게 dose를 정리함 (MCNP의 tally 에 cell 번호 다는 것과 비슷)
[dose ID] [dose 이름] [organ ID list]

* effective dose 정리 관련

RunAction.cc 에 SetEffectiveDose() 에 구현되어 있음. 필요시 수정하면 됨
