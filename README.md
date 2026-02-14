
<!DOCTYPE html>
<html lang="ko">
<body>

    <h1>🚀 Comparison of Overall Survival After Open vs. Robot-Assisted Radical Cystectomy</h1>
    <p><strong>Subtitle: A Target Trial Emulation Approach to Mitigate Immortal Time Bias</strong></p>

    <hr>

    <h2>📝 Overview</h2>
    <p>
        본 프로젝트는 국민건강보험공단(NHIS) 빅데이터를 활용하여 방광암 환자의 <strong>개복 근치적 방광절제술(ORC)</strong>과 <strong>로봇 보조 근치적 방광절제술(RARC)</strong> 간의 전체 생존율(Overall Survival)을 분석한 연구입니다. 
        관찰 연구의 고질적인 문제인 <strong>불멸 시간 편향(Immortal Time Bias)</strong>을 해결하기 위해 <strong>Target Trial Emulation</strong> 프레임워크를 적용하여 인과 추론의 정확도를 높였습니다.
    </p>

    <hr>

    <h2>🔍 Research Background</h2>
    <ul>
        <li><strong>임상적 배경:</strong> 근치적 방광절제술은 근침윤성 방광암의 표준 치료이며, 최근 RARC의 도입이 확대되고 있습니다.</li>
        <li><strong>방법론적 한계:</strong> 기존 관찰 연구들은 수술 시점까지 생존한 환자들만 분석에 포함되는 '불멸 시간 편향'과 환자 상태에 따른 '선택 편향'으로 인해 결과의 왜곡이 발생해 왔습니다.</li>
        <li><strong>연구 목적:</strong> 가상 임상 시험 설계를 통해 이러한 편향들을 보정하고, 두 수술 방식 간의 진정한 생존 혜택 차이를 규명하고자 했습니다.</li>
    </ul>

    <hr>

    <h2>✨ Key Contributions</h2>
    <ul>
        <li><strong>빅데이터 코호트 구축:</strong> 2007~2020년 NHIS 데이터를 활용하여 91,473명의 환자 중 엄격한 기준을 통과한 45,595명의 대규모 코호트를 분석했습니다.</li>
        <li><strong>인과 추론 프레임워크 구현:</strong> 무작위 배정 임상 시험의 논리를 관찰 데이터에 이식하기 위해 적격성, 치료 전략, 타임 제로(Time Zero)를 정교하게 정의했습니다.</li>
        <li><strong>복합 편향 보정:</strong> Cloning 기법으로 불멸 시간 편향을 제거하고, IPCW를 통해 정보성 중도 절단에 따른 선택 편향을 보정했습니다.</li>
    </ul>

    <hr>

    <h2>🛠️ Modeling Strategy & Tech Stack</h2>
    
    <h3>1) Data Handling & Preprocessing (SAS)</h3>
    <ul>
        <li>NHIS 대규모 데이터베이스에서 SQL 기반의 데이터 추출 및 정제 작업을 수행했습니다.</li>
        <li>복잡한 제외 기준(Prior cancer, Missing info, No TUR-BT 등)을 적용하여 분석용 최종 데이터셋을 구축했습니다.</li>
    </ul>

    <h3>2) Causal Inference & Modeling (R)</h3>
    <ul>
        <li><strong>Cloning-Censoring-Weighting (CCW):</strong> 모든 환자를 양측 전략에 동시 배정하는 알고리즘을 구현했습니다.</li>
        <li><strong>IPCW (Inverse Probability of Censoring Weights):</strong> Cox 모델을 활용하여 공변량을 반영한 안정화된 가중치를 산출하고 적용했습니다.</li>
        <li><strong>G-computation:</strong> 표준화(Standardization)를 통해 모집단 수준의 평균 생존 확률 곡선을 시각화했습니다.</li>
        <li><strong>Bootstrap:</strong> 500회 리샘플링을 통해 추정치의 통계적 유의성을 검증했습니다


