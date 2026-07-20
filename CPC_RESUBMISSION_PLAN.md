# CPC 재투고 계획 (LAMMPS-CO)

작성일: 2026-07-20 · 대상 저널: **Computer Physics Communications (CPC)**, Computer Programs in Physics (CPiP)

---

## 1. 투고 이력과 리젝 사유

| 시점 | 저널 | 결과 |
|---|---|---|
| 2026-06-30 | npj Computational Materials | 반려 |
| 2026-07-02 | Archives of Computational Methods in Engineering (ACME) | **Desk reject — scope 불일치** |

ACME 편집장 통지문:

> "The content of the paper is interesting but it falls outside the aims and objectives of the journal,
> focusing on state of the art reviews on computational methods in engineering."

**중요: 동료심사에 가지 않았다.** ACME는 state-of-the-art 리뷰 저널이고 본 원고는 리뷰가 아니므로
편집장 단계에서 반려된 것이다. 따라서 **원고의 학술적 품질에 대한 리뷰어 피드백은 존재하지 않으며,
"어떤 내용이 부족해서 떨어졌다"는 해석은 근거가 없다.**

---

## 2. CPC 필수 요건 (2026-07-20 확인)

CPiP는 "CPC Program Library에 등재될 신규 또는 개정 프로그램을 기술하는" 논문 유형이다.
즉 CPC는 **소프트웨어 저널**이며, 원고보다 **코드 공개와 기탁**이 형식 요건의 중심이다.

| 요건 | 현재 상태 |
|---|---|
| **Program Summary 절** (초록 직후) | ❌ 없음 |
| **Mendeley Data 기탁** → Program Files DOI | ❌ 미기탁 |
| **승인된 오픈소스 라이선스** | ⚠️ LAMMPS가 GPL-2.0 → 승인 목록에 포함되므로 문제 없음 |
| **공개 개발 저장소 링크** | ❌ 원고에 "being prepared for public release" — **투고 전 반드시 공개** |
| 문서화 / User Manual | ⚠️ 미비 |
| 참고문헌 | ❌ "to add on venue selection" 상태 |

승인 라이선스 목록: MIT, Apache-2.0, BSD 3-clause, BSD 2-clause, **GPLv3, GPLv2**, LGPL, MPL-2.0, CeCILL, CeCILL-B.

### Program Summary 필드 구성

CPC 게재 논문 실물(Comput. Phys. Commun. 246 (2020) 106854) 기준:

```
Program Title:
Program Files doi:
Licensing provisions:
Programming language:
Nature of problem:
Solution method:
```

### CPC의 프로그램 채택 기준

> 다른 물리학자에게 유익하거나, 좋은 프로그래밍 관행의 모범이거나,
> **계산물리 커뮤니티에 중요한 새롭거나 참신한 프로그래밍 기법을 보여줄 것.**
> 널리 사용 가능한 언어·하드웨어에서 실행 가능하고, 충분히 문서화되어 있을 것.

본 연구는 세 번째 조항(새로운 프로그래밍 기법)에 정면으로 해당한다.

---

## 3. ⚠️ 최대 위험 — scope 재실패

현재 제목은 **LLM 방법론 논문**의 제목이다.

> "AI-Assisted Systematic Optimization of the LAMMPS Molecular Dynamics Simulator
> **via Large Language Model Coding Agents**"

CPC 편집장이 "계산물리 소프트웨어가 아니라 AI 도구 연구"로 읽을 여지가 있으며,
**ACME와 동일한 유형의 desk reject가 재발할 수 있다.**

### 재프레이밍 방침

| | 현재 | CPC용 |
|---|---|---|
| 주역 | LLM 코딩 에이전트의 능력 | **소프트웨어 산출물** (OMP pair style 35종, 최적화 102파일) |
| LLM | 연구 대상 | 산출물을 만든 **방법** → Methods 절로 이동 |
| 검증 | 부수적 | **핵심 기여** — bit-identical thermo 전수 검증 |
| 응용 사례 | 부록 성격 | **전면 배치** — OPLS-AA TI 63윈도우 ΔG 계산 |

LLM 서사를 버리는 것이 아니다. "전문가 수주 작업을 며칠로 단축"은 여전히 핵심 매력이며
CPC의 "참신한 프로그래밍 기법" 조항에도 부합한다.
다만 **제목과 초록 첫 문장은 소프트웨어여야 한다.**

FEP/TI 사례(21-window × 3용매 = 63윈도우, ΔG_TI = −183.7 / −182.1 / −159.4 kcal/mol)가 특히 중요하다.
CPC가 원하는 "이 소프트웨어로 실제 물리 문제를 풀었다"는 증거가 정확히 이것이다.

---

## 4. 우선순위 (ACME 기준에서 뒤집힘)

| 작업 | ACME 기준 | **CPC 기준** | 근거 |
|---|---|---|---|
| MPI × OpenMP 스케일링 | 2순위 | **1순위** | CPC 독자는 HPC 사용자. 현재 `BUILD_MPI=OFF` |
| Linux/GCC 이식성 검증 | 3순위 | **2순위** | CPC는 "널리 사용 가능한 하드웨어" 요구. LAMMPS 사용자 다수가 Linux/GCC |
| 회귀 테스트 스위트 | 4순위 | **3순위** | 소프트웨어 저널은 검증 가능성을 중시 |
| LLM 재현성 실험 | 1순위 | **선택 사항** | 리젝 사유가 아니었고, CPC에서 방법론은 주역이 아님 |

### 현재 벤치마크의 구조적 약점

전 측정이 **Windows 11 / MSVC 19.43 / AMD Ryzen 9 / `-DBUILD_MPI=OFF`** 단일 환경.
LAMMPS의 본체가 MPI 도메인 분할임을 감안하면, HPC 독자에게는 "실사용 환경 검증 없음"으로 읽힌다.

---

## 5. 실행 계획

### A. 코드 공개 — 계산과 무관, 최우선

저장소가 비공개인 한 CPC 투고 자체가 불가능하다. 전체 일정의 병목일 가능성이 높다.

1. GitHub 저장소 공개 (**GPL-2.0** 명시 — LAMMPS 파생물이므로 필수이자 CPC 승인 목록 포함)
2. 35종 각각에 LAMMPS 스타일 문서(`.rst`) 작성
3. 빌드·재현 절차 문서화 (User Manual 상당)
4. Mendeley Data 기탁 → **Program Files DOI** 발급

### B. iREMB 실험 (계산 비용 가벼움, 하루 안쪽)

32k 원자 500스텝 벤치마크는 한 번에 수 초. 40코어 기준 전체가 하루 안쪽이다.

1. **MPI × OpenMP 하이브리드 스케일링** — 2노드 40코어
   - pure MPI 40랭크, 하이브리드 (2×20 / 4×10 / 8×5), 2노드 강·약 스케일링
   - 신규 35종이 MPI 환경에서 정상 동작함을 동시에 입증
2. **Linux/GCC/Haswell에서 bit-identical 재검증**
   - "MSVC 특유의 아티팩트" 반론 봉쇄
   - `/fp:fast`+AVX2 주장을 `-march=core-avx2 -ffast-math`로 교차검증
3. **회귀 테스트 스위트**
   - 35종 × (1/2/4/8/20 스레드 × MPI 랭크) thermo bit-identity
   - 전 스타일 10k 스텝 이상 NVE 에너지 드리프트 (현재는 2종 × 250스텝뿐)
   - CI 스크립트로 배포 → "자동 수치 검증이 방법론의 성립 조건"이라는 결론을 실증

### C. 원고 재구성

- Program Summary 절 신설 (§2의 필드 구성)
- 제목·초록 재프레이밍 (§3)
- **Limitations 절 신설** — 현재 독립 절이 없고, MPI 부재·단일 OS/CPU·32k 단일 크기가 어디에도 자인되어 있지 않음
- 참고문헌 완성 (LAMMPS, SLEEF, OPLS-AA, TI/BAR 등)
- **35 vs 36종 수치 일관성 점검** — 원고는 35종, 구두 언급은 36종

---

## 6. 선택 사항: LLM 재현성 실험

리젝 사유가 아니었으므로 critical path에서 제외하되, 진짜 약점이므로 여유가 되면 수행한다.

**문제**: "35종 중 27종(77%) 1회 통과"가 **단일 세션 단일 표본**이며,
temperature·sampling 파라미터 미보고, "persistent in-context session"이라 원리적으로 재현 불가한 설정인데
그 사실이 원고에서 논의되지 않는다.

**근거로 쓸 수 있는 공식 문서 사실**: Anthropic 마이그레이션 가이드가 명시하기를,
`temperature = 0`도 이전 모델들에서 동일 출력을 보장한 적이 없다.
따라서 목표는 "재현된다"가 아니라 **"LLM 코드 생성은 비트 단위로 재현되지 않으므로 성공률은 분포로 보고해야 한다"**.

### 설계

| 항목 | 값 |
|---|---|
| 대상 | 35종 전체 |
| 반복 | N=5 독립 세션 (컨텍스트 공유 없음) → 175회 생성 |
| 모델 | `claude-sonnet-4-6` — 원고와 동일 모델 (현재도 서비스 중) |
| 하네스 | Claude Agent SDK / Claude Code CLI — **원고와 동일해야 함** |
| Arm A | 파라미터 미지정 (원고 조건) |
| Arm B | `temperature=0` — 결정론 시도에도 출력이 갈린다는 대조군 |

Sonnet 4.6은 `temperature`를 받는다(Sonnet 5·Opus 4.7 이상에서는 제거되어 400 반환 —
이 사실 자체도 "결정론 제어가 API에서 사라지는 추세"라는 논거가 된다).

### 측정 지표

1. 1회 통과율 — 평균 ± 표준편차, Wilson CI
2. 생성물 간 코드 구조 유사도 (AST/diff)
3. **5개 생성물이 모두 bit-identical thermo를 내는가** ← 가장 중요
4. 생성물별 speedup 편차

3번이 핵심이다. "코드는 매번 다르지만 수치 결과는 항상 같다"가 입증되면
비결정성에도 불구하고 방법론이 신뢰 가능하다는 강력한 논거가 된다.

### 비용

생성당 입력 ~60k · 출력 ~8k 토큰 가정, Sonnet 4.6 단가 $3/$15 per MTok:

| | 캐싱 없이 | 프롬프트 캐싱 적용 |
|---|---|---|
| Arm A | ~$53 | ~$31 |
| A+B 합계 | ~$106 | **~$62** |

같은 style의 5회 반복은 프리픽스가 동일하므로 캐시 쓰기 1회(1.25×) + 읽기 4회(0.1×) = 평균 0.33×.
**주의**: 프리픽스가 바이트 단위로 동일해야 한다. 프롬프트에 타임스탬프·UUID·세션 ID가 들어가면
캐시가 전부 무효화된다. `usage.cache_read_input_tokens`로 검증할 것.
Batch API(50% 할인)는 에이전트 루프에 적용 불가.

### 실행 위치

| 단계 | 위치 | 이유 |
|---|---|---|
| 코드 생성 175회 | **로컬 PC** | LLM API 호출 — iREMB 계산 노드는 외부 네트워크 차단 |
| 빌드 + 50스텝 검증 | 로컬 PC | pair style 1개 재컴파일+링크 ≈ 30초, 총 ~1.5시간 |
| 성능 벤치마크 | **iREMB** | 균질 20코어, 배치 스케줄러로 노이즈 없는 측정 |

---

## 7. 미해결 사항

- **코드 공개 준비 상태** — 라이선스 헤더, 문서, 정리 수준에 따라 A단계 소요가 크게 달라짐
- 백업 저널: **SoftwareX** (같은 Elsevier, 소프트웨어 기술 논문, 심사 부담 적음).
  다만 임팩트와 커뮤니티 적합도는 CPC가 낫다.

---

## 참고

- [Guide for authors — Computer Physics Communications (Elsevier)](https://www.sciencedirect.com/journal/computer-physics-communications/publish/guide-for-authors)
- [CPC International Program Library on Mendeley Data](https://journals.elsevier.com/computer-physics-communications/cpc-international-program-library-on-mendeley-data)
- [예시 CPiP 논문 — AhKin, Comput. Phys. Commun. 246 (2020) 106854](https://mtao8.math.gatech.edu/papers/20CPC.pdf)
