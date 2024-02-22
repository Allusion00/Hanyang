% 두 그룹을 하나의 3D 플롯으로 결합하기

% figure 생성
figure;

% group1을 그래프로 플롯
surf(group1);
hold on; % 다음 플롯을 같은 축에 추가

% group2를 그래프로 플롯
surf(group2);

% 그래프 제목과 범례
title('2D2G Flux Distribution');
legend('Group 1', 'Group 2');

% 그리드 표시
grid on;

% hold 해제
hold off;
