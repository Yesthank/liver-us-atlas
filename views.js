/* 17 standard liver/upper-abdomen ultrasound views.
 * Geometry in BodyParts3D mm (+X patient left, +Y posterior, +Z cranial).
 *   E : approximate skin contact (probe snapped to skin along the beam), or
 *   D : beam direction (probe snapped to skin by marching back from T)
 *   T : point the central beam passes through
 *   H : direction shown on the LEFT of the screen (projected into the scan plane)
 *   normalAxis : for short-axis views, the plane is made perpendicular to this axis
 */
window.US_VIEWS = [
  {
    id: 1, en: 'Substernal longitudinal scan', ko: '검상돌기하 종주사',
    E: [22, -216, 1128], T: [22, -140, 1140], H: [0, 0, 1], depth: 150,
    probe: '검상돌기 바로 아래, 정중선에서 약간 왼쪽. 프로브를 세로로 세우고 마커는 머리 쪽.',
    breath: '숨을 깊게 들이마신 뒤 참기',
    check: ['간 좌외측구역(S2·S3)의 표면과 하연(edge)', '정맥관인대 열 뒤의 미상엽(S1)', '간 뒤쪽을 지나는 대동맥과 복강동맥 분지', '좌엽 크기와 모서리 각도로 간경변 징후 평가'],
  },
  {
    id: 2, en: 'Substernal transverse scan', ko: '검상돌기하 횡주사',
    E: [8, -217, 1132], T: [10, -150, 1160], H: [-1, 0, 0], depth: 150,
    probe: '검상돌기 아래 정중선에 가로로 대고 머리 쪽으로 기울임.',
    breath: '깊은 흡기 후 참기',
    check: ['좌문맥 제대부(umbilical portion)와 분지', '제대열·원인대를 경계로 S3와 S4 구분', '좌간정맥을 따라 S2·S3 구분', '좌엽 실질의 결절 여부'],
  },
  {
    id: 3, en: 'Rt. Subcostal scan ① (confluence)', ko: '우늑골하 주사 ① 간정맥 합류부',
    E: [-35, -214, 1092], T: [-18, -122, 1166], H: [-1, 0, -0.35], depth: 170,
    probe: '오른쪽 늑골연 바로 아래, 늑골연과 나란히 비스듬히 대고 머리 쪽으로 크게 기울임.',
    breath: '깊은 흡기 후 참기',
    check: ['우·중·좌간정맥이 하대정맥으로 합류하는 모습', '우간정맥: S7(뒤)와 S8(앞) 경계', '중간정맥: S8과 S4 경계(Cantlie 선)', '간정맥 구경과 주행(울혈·Budd–Chiari 감별)'],
  },
  {
    id: 4, en: 'Rt. Subcostal scan ② (HV level)', ko: '우늑골하 주사 ② 간정맥 높이',
    E: [-55, -212, 1090], T: [-45, -108, 1146], H: [-1, 0, -0.35], depth: 170,
    probe: '①보다 약간 바깥쪽에서 덜 기울여 간정맥 몸통을 길게 봄.',
    breath: '깊은 흡기 후 참기',
    check: ['간정맥 주행으로 우전구역(S8)과 우후구역(S7) 구분', '중간정맥 좌우의 S4와 S8', '간 실질 에코와 균일도'],
  },
  {
    id: 5, en: 'Rt. Subcostal scan ③ (Dome level)', ko: '우늑골하 주사 ③ 간 돔',
    E: [-72, -209, 1090], T: [-62, -110, 1178], H: [-1, 0, -0.3], depth: 180,
    probe: '오른쪽 늑골연 아래에서 머리 쪽으로 최대한 기울여 횡격막 바로 아래까지.',
    breath: '최대 흡기 후 참기',
    check: ['간 돔(S7·S8 상부)과 횡격막 하부: 놓치기 쉬운 사각지대', '횡격막 위 폐에 의한 거울상 인공물', '우측 흉막 삼출 유무'],
  },
  {
    id: 6, en: 'Rt. Subcostal scan ④ (GB level)', ko: '우늑골하 주사 ④ 담낭 높이',
    E: [-62, -213, 1074], T: [-50, -135, 1062], H: [-1, 0, -0.3], depth: 160,
    probe: '늑골연 아래에서 기울기를 줄여 간 하부와 담낭을 지나게 함.',
    breath: '흡기 후 참기',
    check: ['간 하부 S5·S6', '담낭과 담낭와(GB fossa)', '간 하연과 우신 상극의 관계'],
  },
  {
    id: 7, en: 'Rt. Subcostal scan ⑤ (RPV)', ko: '우늑골하 주사 ⑤ 우문맥',
    E: [-60, -213, 1086], T: [-48, -125, 1120], H: [-1, 0, -0.3], depth: 160,
    probe: '늑골연 아래 비스듬히, 간문부를 향해 중간 정도 기울임.',
    breath: '흡기 후 참기',
    check: ['우문맥이 전분지(S5·S8)와 후분지(S6·S7)로 갈라지는 모습', '문맥면: 위(S7·S8)와 아래(S5·S6) 경계', '미상엽과 하대정맥'],
  },
  {
    id: 8, en: 'Rt. Intercostal scan (anterior)', ko: '우늑간 주사 (전방)',
    T: [-55, -140, 1100], D: [0.85, 0.5, 0.05], H: [0, 0.8, 0.6], depth: 160, avoidRib: true,
    probe: '전액와선 부근 오른쪽 늑간, 늑골과 나란히.',
    breath: '편한 호흡 또는 가벼운 흡기',
    check: ['우문맥 전분지와 S5·S8', '담낭과 간문부', '늑골 음영(acoustic shadow)을 피해 늑간에 정확히 위치'],
  },
  {
    id: 9, en: 'Rt. Intercostal scan (posterior)', ko: '우늑간 주사 (후방)',
    T: [-70, -100, 1112], D: [1, 0.05, 0.05], H: [0, 0.8, 0.6], depth: 160, avoidRib: true,
    probe: '중액와선 부근 오른쪽 늑간, 늑골과 나란히.',
    breath: '편한 호흡 또는 가벼운 흡기',
    check: ['우간정맥과 우문맥 후분지', '우후구역 S6·S7', '비만하거나 간이 위축된 환자에서 우엽을 보는 주된 창'],
  },
  {
    id: 10, en: 'Rt. Flank scan', ko: '우측 옆구리 관상 주사',
    E: [-135, -165, 1040], T: [-62, -90, 1062], H: [0, 0.4, 0.9], depth: 160,
    probe: '전·중액와선 사이 오른쪽 옆구리, 늑골연 바로 아래에 세로(관상 사위면)로 대고 뒤쪽·머리 쪽으로 기울임. 마커는 머리 쪽.',
    breath: '흡기 후 참기',
    check: ['간–우신 경계(Morison 와): 복수 확인', '간 S6·S7과 우신 피질 에코 비교', '횡격막과 늑막각(흉막 삼출)'],
  },
  {
    id: 11, en: 'Lt. intercostal scan (splenic scan)', ko: '좌늑간 주사 (비장)',
    T: [84, -90, 1085], D: [-0.84, -0.52, 0.1], H: [0, 0.75, 0.65], depth: 140, avoidRib: true,
    probe: '왼쪽 후액와선 부근 늑간, 비장 장축을 따라.',
    breath: '흡기 후 참기',
    check: ['비장 장축 길이(비장비대)', '비문부와 비정맥', '좌측 횡격막·좌신 상극'],
  },
  {
    id: 12, en: 'Transverse scan (Body)', ko: '횡주사 (췌장 체부)',
    T: [12, -124, 1052], D: [0, 1, 0.12], H: [-1, 0, 0], depth: 150,
    probe: '상복부 정중선에 가로로, 검상돌기와 배꼽 사이.',
    breath: '흡기로 간을 내려 음향창으로 사용',
    check: ['췌장 체부와 뒤쪽의 비정맥', '상장간막동맥·대동맥·좌신정맥', '주췌관 확장 여부'],
  },
  {
    id: 13, en: 'Oblique scan (Head)', ko: '사위 주사 (췌장 두부)',
    T: [-28, -128, 1033], D: [0, 1, 0.12], H: [-1, 0, -0.25], depth: 150,
    probe: '상복부 정중선 오른쪽, 오른쪽 끝을 약간 내려 비스듬히.',
    breath: '흡기 후 참기',
    check: ['췌장 두부와 십이지장', '상장간막정맥·문맥 합류부', '하대정맥과 원위 총담관'],
  },
  {
    id: 14, en: 'Oblique scan (Tail)', ko: '사위 주사 (췌장 미부)',
    T: [48, -108, 1064], D: [-0.1, 1, 0.15], H: [-0.87, -0.25, -0.42], depth: 160,
    probe: '상복부 왼쪽, 왼쪽 끝을 올려 췌장 장축(비문부 방향)을 따라.',
    breath: '흡기, 필요하면 물을 마셔 위를 음향창으로',
    check: ['췌장 미부와 비정맥', '비문부와 좌신', '위 내 가스에 가리는지 확인'],
  },
  {
    id: 15, en: 'External bile duct longitudinal scan', ko: '간외담관 종주사',
    T: [-24, -132, 1075], D: [0.07, 0.87, 0.33], H: [0.18, -0.36, 0.92], depth: 150,
    probe: '오른쪽 늑골연 아래 정중선 근처, 간십이지장인대를 따라 비스듬히(오른쪽 어깨 방향).',
    breath: '흡기 후 참기, 필요하면 좌측와위',
    check: ['총간관·총담관 구경 측정(대개 6 mm 이하)', '바로 뒤를 지나는 주문맥', '담관과 문맥 사이의 간동맥 단면'],
  },
  {
    id: 16, en: 'GB long axis scan', ko: '담낭 장축 주사',
    T: [-47, -155, 1069], D: [-0.2, 0.75, -0.35], H: [0.31, 0.76, 0.57], depth: 140,
    probe: '오른쪽 늑골연 아래, 담낭 장축(경부→기저부)을 따라.',
    breath: '흡기 후 참기, 좌측와위 추가',
    check: ['경부·체부·기저부 전체', '담석(음영 동반, 이동성)·용종', '벽 두께(3 mm 이하)'],
  },
  {
    id: 17, en: 'GB short axis scan', ko: '담낭 단축 주사',
    T: [-47, -155, 1069], D: [-0.2, 0.75, -0.35], H: [-1, 0, 0], normalAxis: [0.31, 0.76, 0.57], depth: 140,
    probe: '장축 위치에서 프로브를 90° 돌려 담낭을 가로로 자름.',
    breath: '흡기 후 참기',
    check: ['담낭의 둥근 단면', '벽 두께와 벽 주위 액체', '경부에 박힌 결석 확인'],
  },
];
