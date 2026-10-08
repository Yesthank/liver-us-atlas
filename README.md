# 간 초음파 3D 아틀라스

표준 간·상복부 초음파 17개 영상을 세 가지로 함께 보여 주는 교육용 웹페이지입니다.

- 프로브가 몸 표면으로 이동하는 애니메이션
- three.js 3D 해부 모델(Couinaud 분절 S1–S8, 문맥·간정맥·하대정맥·동맥, 담낭·담관, 췌장, 비장, 신장, 소화관, 늑골)과 그 위의 스캔 단면
- 같은 단면으로 합성한 B-mode 초음파 영상과 구조물 라벨

## 여는 법

`index.html`을 더블클릭해 엽니다(인터넷 연결 필요: three.js와 글꼴을 CDN에서 불러옵니다).
브라우저가 로컬 파일 접근을 막으면 `serve.bat`을 실행하세요. 로컬 서버(포트 8917)가 켜지고 브라우저에서 `http://localhost:8917/`가 열립니다.

## 조작

| 동작 | 방법 |
|---|---|
| 영상 선택 | 왼쪽 목록 클릭, 또는 ← → |
| 3D 회전·확대·이동 | 드래그 · 휠 · 우클릭 드래그 |
| 프로브 다시 대기 | `다시 대기` |
| 단면 훑기 | `스윕` 또는 Space |
| 17개 자동 재생 | `순서대로 재생` |
| 미세 조정 | 기울기(fanning) · 회전 · 밀기 슬라이더 |
| 초음파 설정 | 구조 이름, 분절 색 덧씌우기, 게인, 깊이 |

## 실제 초음파 영상 넣기

`us_images/01.jpg` ~ `us_images/17.jpg`(또는 `.png`)를 넣으면 해당 영상 화면의 "실제 초음파 영상" 칸에 자동으로 나타납니다.
그 칸에 이미지 파일을 끌어다 놓아도 되며, 이 경우 그 브라우저에만 저장됩니다.

## 파일

```
index.html          앱 본체 (3D, 초음파 합성, UI)
views.js            17개 표준영상 정의: 프로브 위치·빔 방향·화면 방향, 검사 방법, 확인할 것
data/anatomy.js     압축된 메시 40개 + 1.5 mm 라벨 볼륨 (BodyParts3D 파생 데이터)
tools/              데이터 재생성 스크립트
serve.bat           로컬 서버 실행
```

영상 위치를 고치려면 `views.js`만 수정합니다. 각 영상은 다음 값으로 정의합니다.
`T`는 중심 빔이 지나는 점, `E`는 피부 접촉 근사점이고 `E`가 없으면 `D`(빔 방향)를 씁니다. `H`는 화면 왼쪽에 올 방향입니다. 좌표는 BodyParts3D mm(+x 환자 좌측, +y 등쪽, +z 머리쪽)입니다.

## 데이터 다시 만들기

Python 3.11, `numpy scipy trimesh fast-simplification scikit-image`가 필요합니다. 작업 폴더에 `groups.json`을 복사한 뒤 차례로 실행합니다.

```bash
python tools/fetch_index.py      # BodyParts3D isa zip 목차만 HTTP Range로 읽기 (zip_index.json)
python tools/fetch_objs.py       # groups.json에 적힌 OBJ만 받아 obj/ 에 저장 (~15 MB)
python tools/build_volume.py obj .
python tools/build_assets.py . data/anatomy.js
python tools/portal_segments.py data/anatomy.js   # 구역은 간정맥·열, 분절은 문맥 영역으로 다시 나눔
```

## 출처와 한계

- 3D 해부: BodyParts3D, © The Database Center for Life Science, CC BY-SA 2.1 JP. `data/anatomy.js`도 같은 라이선스를 따릅니다.
- Couinaud 분절: 원본 분절 메시는 라벨이 실제 위치와 어긋나 있어 `tools/portal_segments.py`로 다시 계산했습니다. Couinaud처럼 구역(sector)은 간정맥이 지나는 면으로, 분절은 문맥 분지가 맡는 영역으로 나눴고, 이름은 Brisbane 2000 용어를 따릅니다.
  - 좌우 간(Cantlie 선): 담낭와에서 합류부 직전의 중간정맥 본간까지 잇는 면입니다.
  - 우전(S5·S8)·우후(S6·S7) 구역: 높이마다 우간정맥이 지나는 면입니다. S7·S8 문맥 분지 옆을 따라가는 지류는 빼고, 면이 전·후 문맥 분지 사이를 벗어나지 않게 했습니다. 우간정맥이 끝난 위아래 높이에서는 끝의 각도를 그대로 이어 씁니다.
  - S2·S3: 좌간정맥이 지나는 면(Couinaud의 좌측 문맥 열)으로 나눕니다.
  - S3·S4: 좌문맥 제대부를 지나는 제대열(겸상인대·원인대 선)로, 하대정맥 바로 앞에서 약간 왼쪽에 놓입니다.
  - S5/S8, S6/S7, S4a/S4b: 구역 안에서 가장 가까운 분절 문맥 가지의 영역입니다. 분지 위치는 스크립트의 `SEEDS`에 있습니다.
  - 문맥 이름은 공급하는 쪽으로 붙였습니다(S5–S8은 우문맥, S2–S4는 좌문맥).
  - S1: 원본 미상엽 메시를 사용했습니다.
  - 한계: 이 모델의 문맥 트리는 몇 단계에서 끝나므로 분절 경계는 실제 문맥 영역의 근사입니다. 모델의 혈관 배치 때문에 S5 분지의 일부는 담낭와 높이에서 S4b 쪽에 걸립니다. S2·S3 문맥 분지는 둘 다 좌간정맥보다 앞쪽을 지나므로, 좌간정맥 면으로 나누면 S2 분지의 대부분이 S3 쪽에 놓입니다(문맥 영역으로 나누면 반대로 좌간정맥이 S2 안에 묻힙니다). S8 앞쪽 가지 하나는 복셀에서 좌문맥과 맞닿아 있어 골격 트리상 좌문맥에 붙고, 그 앞부분은 S4 안을 지납니다. S7은 문맥 가지가 없는 후하부까지 내려옵니다. 평면 기준 분절과 실제 문맥 영역이 넓게 어긋날 수 있다는 점도 알려져 있습니다(Fasel 외 1998). 실제 환자의 변이와 다를 수 있습니다.
  - 참고: Couinaud C. *Le Foie: études anatomiques et chirurgicales*. Masson, 1957 · Bismuth H. World J Surg 1982;6:3–9 · Strasberg SM 외. The Brisbane 2000 Terminology of Liver Anatomy and Resections. HPB 2000;2:333–339 · Soler L 외. Comput Aided Surg 2001;6:131–142 · Selle D 외. IEEE Trans Med Imaging 2002;21:1344–1357 · Fasel JH 외. Radiology 1998;206:151–156.
- 초음파 영상: 조직별 에코, 감쇠, 경계면 반사, 스펙클, 뼈·가스 음영, 후방 증강, 횡격막 거울상을 단순 규칙으로 흉내 낸 시뮬레이션입니다.
- 표준영상 구성: 사용자가 정한 17개 목록을 따랐습니다. 참고 자료는 대한초음파의학회 상복부 초음파 표준영상(일반 12·정밀 15)과 국립암센터 간암 검진 질지침입니다.
