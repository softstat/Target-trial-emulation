
# 개복 vs 로봇 보조 방광절제술 후 전체 생존기간 비교  
## Target Trial Emulation을 활용한 Immortal Time Bias 완화

> **Journal**: The Korean Journal of Applied Statistics (KJAS)  
> **Manuscript ID**: KJAS-2025-052  
> **Authors**: 김민석, 이주영 (중앙대학교 응용통계학과)  
> **사용 언어**: SAS, R

---

## 1. 연구 배경 및 목적

### 배경
- **방광암(Bladder Cancer)**은 전 세계적으로 흔한 악성 종양이며, **근치적 방광절제술(Radical Cystectomy)**이 근침윤성 방광암의 표준 치료법
- 2000년대 초반부터 **로봇 보조 근치적 방광절제술(RARC)**이 최소 침습 대안으로 부상
- 기존 관찰 연구들은 **Immortal Time Bias**와 **Selection Bias**로 인해 일관되지 않은 결과를 보고

### 목적
- **Target Trial Emulation** 프레임워크를 적용하여 **ORC(개복)** vs **RARC(로봇)** 간 **전체 생존율(OS)** 비교
- 기존 분석 방법(Naive 분석)과 **Cloning-Censoring-Weighting(CCW)** 방법의 결과를 비교하여 편향의 영향을 검증

---

## 2. 핵심 개념 정리

### 2.1 Immortal Time Bias (불멸 시간 편향)
| 항목 | 설명 |
|------|------|
| **정의** | 적격 기준 시점(TUR-BT)과 치료 배정 시점(수술) 사이에 사망이 발생할 수 없는 기간이 존재하여 생존 추정치가 왜곡되는 편향 |
| **발생 원인** | 추적 관찰이 진단/TUR-BT에서 시작되지만, 치료 그룹은 이후 수술 유형으로 분류 → 수술까지 생존해야 그룹에 포함됨 |
| **영향** | 치료 전략에 인위적인 생존 이점 부여 → Hazard Ratio 추정치 왜곡 |

### 2.2 Target Trial Emulation
- 관찰 데이터에서 **가상의 무작위 배정 임상시험을 모사**하여 인과 추론을 수행하는 프레임워크
- **CCW(Cloning-Censoring-Weighting)** 방법이 핵심 구현 기법

### 2.3 CCW 방법의 3단계

```
[Cloning]  → 각 환자를 2개 복제: ORC clone + RARC clone
[Censoring] → 실제 치료가 배정 전략에서 이탈하면 해당 clone을 인위적 중도절단
[Weighting] → 인위적 중도절단에 의한 편향을 IPCW로 보정
```

---

## 3. 분석 방법론 비교

### 3.1 전통적 분석 방법 (Naive Approaches)

| 방법 | Time Zero | 치료 배정 | 문제점 |
|------|-----------|---------|--------|
| **Immortal-time-included** | TUR-BT | 미래 수술 정보 사용 | Immortal Time Bias 직접 포함 |
| **Immortal-time-excluded** | 수술일 | 수술 시점 기준 | Selection Bias (수술 전 사망자 제외) |
| **Landmark Analysis** | 2년 시점 | Landmark 전 수술 기준 | Selection Bias (Landmark 전 사망/중도절단자 제외) |

### 3.2 CCW 방법 (Target Trial Emulation)

| 항목 | 내용 |
|------|------|
| **Time Zero** | TUR-BT 시점 |
| **치료 배정** | Cloning을 통해 동시에 양 전략에 배정 (미래 정보 불사용) |
| **Grace Period** | TUR-BT 후 2년 이내 수술 시행 |
| **편향 보정** | IPCW (Inverse Probability of Censoring Weights) |
| **신뢰구간** | Nonparametric Bootstrap (500회) |
| **생존 확률 추정** | G-computation |

### 3.3 통계 모형

**Cox Proportional Hazards Model** (공통 기반):

```
h(t | Z_i) = h_0(t) × exp(β'Z_i)
```

**IPCW 가중 부분우도 (CCW 전용)**:

```
L_IPCW(β) = ∏ [ sw_i(t) × exp(β_A × A_i^g) / Σ sw_j(t) × Y_j(t) × exp(β_A × A_j^g) ]^ΔN_i(t)
```

- `sw_i(t)`: 안정화 역확률 중도절단 가중치 = Ĝ₀(t) / Ĝ(t|z_i)
- 중도절단 확률은 Cox 모형으로 추정 (공변량: 연령, 성별, 동반질환, 진단연도)

---

## 4. 데이터

### 4.1 데이터 출처
- **국민건강보험공단(NHIS) 데이터베이스**: 약 5천만 명의 한국 인구를 포괄
- **연구 기간**: 2007–2020년 (최대 14년 추적 관찰)

### 4.2 대상자 선정 흐름

```
방광암 신규 진단 환자: 91,473명
  ↓ 제외: 타 암 병력(45,350), 결측(167), TUR-BT 미시행(361)
최종 연구 대상: 45,595명
  ↓ Cloning
  ORC군: 4,877명  |  RARC군: 689명
```

### 4.3 주요 변수

| 구분 | 변수 |
|------|------|
| **노출(Exposure)** | 수술 유형 (ORC vs RARC) |
| **결과(Outcome)** | 전체 사망 (All-cause mortality) |
| **공변량(Covariates)** | 성별, 연령, 진단연도, 소득수준, 고혈압, CCI(찰슨동반질환지수) |

### 4.4 기저 특성 요약 (Table 2)

| 특성 | ORC (N=4,877) | RARC (N=689) | P value |
|------|:---:|:---:|:---:|
| 연령 65-75세 | 38.73% | 40.78% | 0.239 |
| 남성 | 85.11% | 84.76% | 0.852 |
| CCI ≥ 4 | 64.00% | **71.00%** | **<0.001** |
| 고혈압 有 | 51.26% | **69.67%** | **<0.001** |

→ RARC군이 동반질환 부담이 유의하게 높음

---

## 5. 주요 결과

### 5.1 CCW 분석 결과 (Primary)

| 그룹 | ORC 생존률 | RARC 생존률 | Risk Difference | Hazard Ratio (95% CI) |
|------|:---:|:---:|:---:|:---:|
| **전체** | 72.62% | 81.38% | 8.76% | **1.02 (0.94, 1.11)** |
| 남성 | 76.28% | 78.38% | 2.10% | 1.01 (0.93, 1.12) |
| 여성 | 74.24% | 82.43% | 8.20% | 1.09 (0.91, 1.30) |

**→ HR ≈ 1.02로 두 수술법 간 유의한 생존 차이 없음**

### 5.2 전통적 분석 결과 (Conventional)

| 방법 | 전체 HR (95% CI) | 해석 |
|------|:---:|------|
| Immortal-time-included | 1.04 (0.92, 1.17) | Immortal time이 HR을 **감쇠(attenuate)** |
| Immortal-time-excluded | **1.42 (1.19, 1.71)** | Selection Bias로 HR **과대추정** |
| Landmark | **1.29 (1.08, 1.55)** | Selection Bias로 HR **과대추정** |

### 5.3 결과 해석
- 전통적 방법에서 ORC의 생존 이점으로 보이던 결과는 **방법론적 편향**에 기인
- CCW를 적용하면 편향이 보정되어 두 수술법 간 **차이 소실**
- TUR-BT↔수술 간 평균 간격: ORC 19.6±9.7개월, RARC 21.9±25.8개월 → RARC에서 immortal time이 더 길고 변동 큼

---

## 6. SAS / R 구현 포인트

### R에서의 주요 패키지 및 구현

```r
# 생존 분석 기본
library(survival)    # Cox PH model, Kaplan-Meier
library(survminer)   # 생존 곡선 시각화

# CCW 구현 핵심 단계
# 1. Cloning: 각 환자를 2행으로 복제 (ORC/RARC 배정)
# 2. Censoring: 배정 전략 이탈 시 인위적 중도절단 적용
# 3. IPCW 가중치 계산
#    - 중도절단 모형: coxph(Surv(time, censor_event) ~ age + sex + cci + dx_year)
#    - 안정화 가중치: KM marginal / Cox conditional
# 4. 가중 Cox 모형: coxph(..., weights = sw)
# 5. G-computation: 조건부 생존함수의 경험 분포 평균

# Bootstrap (500회)
boot_results <- replicate(500, {
  boot_idx <- sample(1:n, replace = TRUE)
  # clone → censor → weight → fit 반복
})
```

### SAS에서의 주요 구현

```sas
/* 생존 분석 */
PROC PHREG DATA=cloned_data;
  CLASS treatment age_group sex cci dx_year;
  MODEL time*event(0) = treatment age_group sex cci dx_year;
  WEIGHT sw;  /* 안정화 IPCW 가중치 */
RUN;

/* 중도절단 모형 (IPCW 산출용) */
PROC PHREG DATA=cloned_data;
  CLASS age_group sex cci dx_year;
  MODEL time*artificial_censor(0) = age_group sex cci dx_year;
  OUTPUT OUT=censor_prob SURVIVAL=G_hat;
RUN;

/* 기저특성 비교 */
PROC FREQ DATA=study_pop;
  TABLES (age_group sex dx_year income cci hypertension hospital)*surgery_type
         / CHISQ;
RUN;
```

---

## 7. 연구 의의 및 한계

### 의의
- **Target Trial Emulation**을 방광암 수술 비교에 적용한 국내 최초 연구
- Immortal Time Bias가 관찰 연구 결과를 어떻게 왜곡하는지 실증적으로 시연
- CCW 방법이 기존 분석 대비 편향을 효과적으로 보정함을 확인

### 한계
- NHIS 데이터에 **종양 병기(T stage), 수행 상태(ECOG), 검사 소견** 등 임상 정보 부재 → 잔여 교란(Residual Confounding) 가능
- RARC 코드 부재로 마취 코드 기반 간접 정의 → 오분류 가능성
- CCW는 **측정된 교란 변수의 완전성**을 가정

---

## 8. 핵심 키워드 정리

| 키워드 | 설명 |
|--------|------|
| **Immortal Time Bias** | 적격 시점~치료 시점 사이 사망 불가능 기간으로 인한 편향 |
| **Target Trial Emulation** | 관찰 데이터에서 RCT를 모사하는 인과 추론 프레임워크 |
| **CCW** | Cloning-Censoring-Weighting; TTE의 핵심 구현 방법 |
| **IPCW** | Inverse Probability of Censoring Weights; 정보적 중도절단 보정 |
| **G-computation** | 조건부 생존함수를 공변량 분포 위에서 표준화하여 인과적 생존 확률 추정 |
| **Grace Period** | 치료 배정이 이루어질 것으로 기대되는 기간 (본 연구: 2년) |
| **Per-protocol Effect** | 배정된 치료를 준수했을 때의 인과 효과 |
| **Cox PH Model** | 공변량의 위험 비례 가정 하에 시간-사건 데이터를 분석하는 준모수 모형 |
| **Bootstrap** | 비모수 재표본 방법; CCW에서 신뢰구간 산출에 활용 |

---

## 9. 면접 대비 Q&A

**Q1. Immortal Time Bias가 왜 발생하나요?**  
> 추적 시작(TUR-BT)과 치료 배정(수술) 시점이 일치하지 않아, 수술까지 생존해야만 해당 치료군에 포함되므로 인위적 생존 이점이 발생합니다.

**Q2. CCW 방법의 핵심 아이디어는?**  
> 모든 환자를 두 치료 전략에 동시에 배정(Cloning)하고, 실제 치료가 배정 전략에서 이탈하면 중도절단(Censoring)한 뒤, 이로 인한 정보적 중도절단을 IPCW로 보정(Weighting)합니다.

**Q3. 왜 안정화 가중치(Stabilized Weights)를 사용하나요?**  
> 역확률 가중치가 극단적으로 커질 수 있어 분산이 증가합니다. 분자에 주변(marginal) 중도절단 생존함수를 넣어 가중치를 안정화하면 분산이 줄어듭니다.

**Q4. Landmark Analysis와 CCW의 차이는?**  
> Landmark은 특정 시점까지 생존한 환자만 포함하므로 Selection Bias가 존재합니다. CCW는 모든 적격 환자를 포함하고 IPCW로 보정하여 TUR-BT 시점에서의 인과 효과를 추정합니다.

**Q5. 이 연구의 주요 결론은?**  
> CCW를 적용하면 ORC와 RARC 간 전체 생존에 유의한 차이가 없었으며(HR=1.02), 기존 연구에서 보고된 ORC의 생존 이점은 방법론적 편향에 기인한 것으로 판단됩니다.


