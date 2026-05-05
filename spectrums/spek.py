import numpy as np
import spekpy
import pandas as pd
from scipy.interpolate import interp1d

# 공통 에너지 축 (0 ~ 100 keV, 0.5 keV 간격)
common_energy = np.arange(0, 101, 0.5)

# 중심 kVp 값들
kvp_centers = [70, 80, 90]
# 필터와 리플 조건: (필터 두께 (mm Al), 리플(%; 0.15는 15%))
filter_conditions = [
    (0.5, 0.15),
    (1.5, 0.15),
    (2.5, 0.10),
    (3.0, 0.05)
]

# 각 조합에 대해 스펙트럼 계산
for mean_kVp in kvp_centers:
    for (filter_thickness, ripple) in filter_conditions:
        # 리플에 따른 전체 전압 변동: 예) 15%의 경우 70 ± (70*0.15)= ±10.5 kVp
        delta = mean_kVp * ripple
        # 5개의 kVp 샘플 (선형 보간)
        sample_kvps = np.linspace(mean_kVp - delta, mean_kVp + delta, 5)
        # 가중치 계산: Gaussian 분포, sigma = delta/2
        sigma = delta / 2.0
        weights = np.exp(-0.5 * ((sample_kvps - mean_kVp) / sigma) ** 2)
        weights /= np.sum(weights)  # 정규화

        spectra_list = []
        
        # 각 kVp 샘플에 대해 스펙트럼 계산 및 보간
        for kvp in sample_kvps:
            s = spekpy.Spek()
            s.set(kvp=kvp)
            s.filter('Al', filter_thickness)
            spectrum = s.get_spectrum()  # 반환: (energy_array, fluence_array)
            energy = spectrum[0]
            fluence = spectrum[1]
            
            # 공통 에너지 축으로 보간 (에너지 bin 사이즈 0.5 keV 적용)
            interp_func = interp1d(energy, fluence, kind='linear', bounds_error=False, fill_value=0)
            interp_fluence = interp_func(common_energy)
            spectra_list.append(interp_fluence)
        
        # 가중 평균 스펙트럼 계산
        final_spectrum = np.average(spectra_list, axis=0, weights=weights)
        # 평균 에너지 계산: Σ(E*fluence) / Σ(fluence)
        avg_energy = np.sum(common_energy * final_spectrum) / np.sum(final_spectrum)
        
        # 결과 출력
        print(f"Mean kVp: {mean_kVp} kVp, Filter: {filter_thickness} mm Al, Ripple: {ripple*100:.0f}% -> Average Energy: {avg_energy:.2f} keV")
        
        # 파일 이름 예: Xray_Spectrum_70kVp_0.5mmAl_15percentRipple.csv
        filename = f"Xray_Spectrum_{mean_kVp}kVp_{filter_thickness}mmAl_{int(ripple*100)}percentRipple.csv"
        df = pd.DataFrame({'Energy (keV)': common_energy, 'Fluence (normalized)': final_spectrum})
        df.to_csv(filename, index=False)

