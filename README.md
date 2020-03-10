# TetCal
(VoxelVersion과 최대한 동기화)

* __아직 구현 안된 부분__: surface source
* __구현된 부분__: External, Internal, Hadronic, dose 정리, effective dose 정리
* __dose 정리 관련__: [팬텀이름].dose 라는 파일 내에 아래 형식에 맞게 dose를 정리함 (MCNP의 tally 에 cell 번호 다는 것과 비슷)
```
[dose ID] [dose 이름] [organ ID list]
```
* __effective dose 정리 관련__: RunAction.cc 에 SetEffectiveDose() 에 구현되어 있음. 필요시 수정하면 됨
