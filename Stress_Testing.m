clear; clc; close all;
D_max = 0.180;
D_with_spikes = 0.184;
D_rotating_region = 0.190;
D_min = 0.158;
D_hub = 0.039;
R_max = D_max / 2;
R_with_spikes = D_with_spikes / 2;
R_rotating = D_rotating_region / 2;
R_min = D_min / 2;
R_hub = D_hub / 2;
omega_rpm = 350;
omega = omega_rpm * 2*pi/60;
v_water = 1.0;
deg_per_step = 3;
dt_cfd = 0.00143;

tip_speed = omega * R_max;
chord_typical = 0.030;
nu_water = 1e-6;
Re_tip = tip_speed * chord_typical / nu_water;

immersion_depth = 0.06;
water_level = -immersion_depth;

fprintf('===== SPECIFICATIONS =====\n');
fprintf('Wheel Max Diameter: %.1f mm\n', D_max*1000);
fprintf('Hub Diameter: %.1f mm\n', D_hub*1000);
fprintf('Rotation: %d RPM (%.2f rad/s)\n', omega_rpm, omega);
fprintf('Tip Speed: %.3f m/s\n', tip_speed);
fprintf('\n');
n_spokes = 6;
hydrofoil_type = 'NACA 0015';
thickness_ratio = 0.15;
n_sections = 30;
r_sections = linspace(R_hub, R_max, n_sections);
chord_hub = 0.040;
chord_tip = 0.025;
chord = linspace(chord_hub, chord_tip, n_sections);

taper_length = R_max - R_hub;
taper_height = (R_max - R_min);
taper_angle = atan(taper_height / taper_length);
span = 2 * (R_min + (r_sections - R_hub) * tan(taper_angle));

pitch_hub = 8;
pitch_tip = 3;
pitch_angle = linspace(pitch_hub, pitch_tip, n_sections);

fprintf('===== HYDROFOIL DESIGN =====\n');
fprintf('Profile: %s\n', hydrofoil_type);
fprintf('Spokes: %d\n', n_spokes);
fprintf('Sections per spoke: %d\n', n_sections);
fprintf('\n');
rho_material = 1140;
E = 3.5e9;
sigma_yield = 85e6;
material_name = 'CF-Nylon';
safety_factor_target = 2.5;
rho_water = 1000;
nu = 1e-6;
g = 9.81;
t_end = 1.0;
dt = dt_cfd;
t = 0:dt:t_end;
n_steps = length(t);

fprintf('Running simulation...\n');
mass_per_spoke = 0;
for i = 1:n_sections-1
    dr = r_sections(i+1) - r_sections(i);
    avg_chord = (chord(i) + chord(i+1)) / 2;
    avg_span = (span(i) + span(i+1)) / 2;
    thickness_local = avg_chord * thickness_ratio;
    volume_section = avg_chord * avg_span * thickness_local * dr / avg_chord;
    mass_per_spoke = mass_per_spoke + rho_material * volume_section;
end
total_spoke_mass = mass_per_spoke;
total_wheel_mass = n_spokes * total_spoke_mass + 0.080;

theta = zeros(n_steps, n_spokes);
F_thrust = zeros(n_steps, n_spokes, n_sections);
F_drag = zeros(n_steps, n_spokes, n_sections);
F_lift_vertical = zeros(n_steps, n_spokes, n_sections);
F_total = zeros(n_steps, 3);
immersed = zeros(n_steps, n_spokes, n_sections);

F_centrifugal = zeros(n_spokes, n_sections);
stress_centrifugal = zeros(n_spokes, n_sections);
stress_bending = zeros(n_spokes, n_sections);
stress_total = zeros(n_spokes, n_sections);
safety_factor = zeros(n_spokes, n_sections);

for i = 1:n_steps
    for j = 1:n_spokes
        theta(i,j) = omega * t(i) + (j-1) * 2*pi/n_spokes;
        
        for k = 1:n_sections
            y_position = r_sections(k) * sin(theta(i,j));
            
            if y_position < water_level
                immersed(i,j,k) = 1;
                depth_below_surface = water_level - y_position;
                immersion_factor = min(1.0, depth_below_surface / (chord(k)));
                
                v_tangential = omega * r_sections(k);
                v_rel_x = -v_tangential * sin(theta(i,j)) + v_water;
                v_rel_y = v_tangential * cos(theta(i,j));
                v_rel_mag = sqrt(v_rel_x^2 + v_rel_y^2);
                
                flow_angle = atan2(v_rel_y, v_rel_x);
                spoke_angle = theta(i,j) + pi/2 + deg2rad(pitch_angle(k));
                alpha = rad2deg(spoke_angle - flow_angle);
                
                while alpha > 180, alpha = alpha - 360; end
                while alpha < -180, alpha = alpha + 360; end
                
                alpha_rad = deg2rad(alpha);
                
                if abs(alpha) < 10
                    Cl = 2*pi * alpha_rad * 0.9;
                    Cd = 0.01 + 0.04 * alpha_rad^2;
                elseif abs(alpha) < 15
                    Cl = 1.3 * sign(alpha);
                    Cd = 0.05 + 0.15 * abs(alpha_rad);
                else
                    Cl = 0.9 * sign(alpha);
                    Cd = 0.4 + 0.4 * abs(sin(alpha_rad));
                end
                
                Re_local = v_rel_mag * chord(k) / nu;
                if Re_local < 1e5
                    Cl = Cl * (0.8 + 0.2 * (Re_local / 1e5));
                end
                
                q = 0.5 * rho_water * v_rel_mag^2;
                
                if k < n_sections
                    dr = r_sections(k+1) - r_sections(k);
                else
                    dr = r_sections(k) - r_sections(k-1);
                end
                A_section = chord(k) * dr;
                
                F_lift = q * A_section * Cl * immersion_factor;
                F_drag_local = q * A_section * Cd * immersion_factor;
                
                F_thrust(i,j,k) = F_lift * cos(flow_angle) - F_drag_local * sin(flow_angle);
                F_drag(i,j,k) = F_lift * sin(flow_angle) + F_drag_local * cos(flow_angle);
                
                if theta(i,j) > pi && theta(i,j) < 2*pi
                    lift_angle = spoke_angle;
                    F_lift_vertical(i,j,k) = F_lift * abs(sin(lift_angle)) * sign(cos(theta(i,j)));
                else
                    F_lift_vertical(i,j,k) = 0;
                end
            else
                immersed(i,j,k) = 0;
                F_thrust(i,j,k) = 0;
                F_drag(i,j,k) = 0;
                F_lift_vertical(i,j,k) = 0;
            end
        end
    end
    
    F_total(i,1) = sum(sum(F_thrust(i,:,:)));
    F_total(i,2) = sum(sum(F_drag(i,:,:)));
    F_total(i,3) = sum(sum(F_lift_vertical(i,:,:)));
end

%% STRESS ANALYSIS
fprintf('Calculating stress distribution...\n');

for j = 1:n_spokes
    for k = 1:n_sections
        if k < n_sections
            dr = r_sections(k+1) - r_sections(k);
        else
            dr = r_sections(k) - r_sections(k-1);
        end
        
        thickness_local = chord(k) * thickness_ratio;
        volume_local = chord(k) * span(k) * thickness_local * dr / chord(k);
        mass_local = rho_material * volume_local;
        
        F_centrifugal(j,k) = mass_local * omega^2 * r_sections(k);
        A_cross = thickness_local * span(k);
        stress_centrifugal(j,k) = F_centrifugal(j,k) / A_cross;
        
        F_hydro_avg = mean(sqrt(F_thrust(:,j,k).^2 + F_drag(:,j,k).^2));
        moment_arm = r_sections(k) - R_hub;
        M_bending = F_hydro_avg * moment_arm;
        
        I = (span(k) * thickness_local^3) / 12;
        c = thickness_local / 2;
        stress_bending(j,k) = M_bending * c / I;
        
        stress_total(j,k) = stress_centrifugal(j,k) + stress_bending(j,k);
        safety_factor(j,k) = sigma_yield / stress_total(j,k);
    end
end

%% RESULTS
fprintf('\n═══════════════════════════════════════\n');
fprintf('Max Centrifugal: %.2f MPa\n', max(stress_centrifugal(:))/1e6);
fprintf('Max Bending:     %.2f MPa\n', max(stress_bending(:))/1e6);
fprintf('Max Total:       %.2f MPa\n', max(stress_total(:))/1e6);
fprintf('Yield:           %.2f MPa\n', sigma_yield/1e6);
fprintf('Min SF:          %.2f\n', min(safety_factor(:)));
fprintf('═══════════════════════════════════════\n\n');

%% HELPER FUNCTION TO DRAW CIRCLES
function draw_circle(center_x, center_y, radius, varargin)
    theta = linspace(0, 2*pi, 100);
    x = center_x + radius * cos(theta);
    y = center_y + radius * sin(theta);
    plot(x, y, varargin{:});
end

%% VISUALIZATION
fprintf('Generating stress visualization...\n');

stress_colormap = jet(256);

% Figure 1: 3D Stress Distribution
figure('Name', 'Total Stress Distribution on Blades', 'NumberTitle', 'off', 'Position', [100 100 1200 800]);

for j = 1:n_spokes
    angle_base = (j-1) * 2*pi/n_spokes;
    
    for k = 1:n_sections-1
        r1 = r_sections(k);
        r2 = r_sections(k+1);
        
        x1 = r1 * cos(angle_base);
        y1 = r1 * sin(angle_base);
        x2 = r2 * cos(angle_base);
        y2 = r2 * sin(angle_base);
        
        half_span1 = span(k) / 2;
        half_span2 = span(k+1) / 2;
        
        X = [x1 x2 x2 x1] * 1000;
        Y = [y1 y2 y2 y1] * 1000;
        Z = [-half_span1 -half_span2 half_span2 half_span1] * 1000;
        
        stress_val = stress_total(j,k);
        stress_normalized = (stress_val - min(stress_total(:))) / (max(stress_total(:)) - min(stress_total(:)));
        color_idx = max(1, min(256, round(stress_normalized * 255 + 1)));
        face_color = stress_colormap(color_idx, :);
        
        patch(X, Y, Z, face_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
end

xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('Total Stress Distribution on Blades (3D View)');
colormap(jet);
c = colorbar;
c.Label.String = 'Stress (MPa)';
caxis([min(stress_total(:))/1e6, max(stress_total(:))/1e6]);
axis equal;
grid on;
view(45, 30);
lighting gouraud;
light('Position', [1 1 1]);

% Figure 2: Top View - Total Stress
figure('Name', 'Top View - Stress Distribution', 'NumberTitle', 'off', 'Position', [100 100 1000 900]);

for j = 1:n_spokes
    angle_base = (j-1) * 2*pi/n_spokes;
    
    for k = 1:n_sections-1
        r1 = r_sections(k);
        r2 = r_sections(k+1);
        
        x1 = r1 * cos(angle_base);
        y1 = r1 * sin(angle_base);
        x2 = r2 * cos(angle_base);
        y2 = r2 * sin(angle_base);
        
        width = chord(k);
        dx = width * sin(angle_base);
        dy = -width * cos(angle_base);
        
        X = [x1, x2, x2+dx, x1+dx] * 1000;
        Y = [y1, y2, y2+dy, y1+dy] * 1000;
        
        stress_val = stress_total(j,k);
        stress_normalized = (stress_val - min(stress_total(:))) / (max(stress_total(:)) - min(stress_total(:)));
        color_idx = max(1, min(256, round(stress_normalized * 255 + 1)));
        face_color = stress_colormap(color_idx, :);
        
        patch(X, Y, face_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
end

draw_circle(0, 0, R_hub*1000, 'k', 'LineWidth', 2);
draw_circle(0, 0, R_max*1000, 'k--', 'LineWidth', 1);

xlabel('X (mm)');
ylabel('Y (mm)');
title('Top View - Total Stress Distribution');
colormap(jet);
c = colorbar;
c.Label.String = 'Total Stress (MPa)';
caxis([min(stress_total(:))/1e6, max(stress_total(:))/1e6]);
axis equal;
grid on;

text_x = R_max * 1.3 * 1000;
text_y = R_max * 1.1 * 1000;
text(text_x, text_y, sprintf('Max: %.1f MPa', max(stress_total(:))/1e6), 'FontSize', 10, 'FontWeight', 'bold');
text(text_x, text_y-10, sprintf('Min: %.1f MPa', min(stress_total(:))/1e6), 'FontSize', 10, 'FontWeight', 'bold');
text(text_x, text_y-20, sprintf('Yield: %.1f MPa', sigma_yield/1e6), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');

% Figure 3: Centrifugal Stress
figure('Name', 'Centrifugal Stress Distribution', 'NumberTitle', 'off', 'Position', [150 150 1000 900]);

for j = 1:n_spokes
    angle_base = (j-1) * 2*pi/n_spokes;
    
    for k = 1:n_sections-1
        r1 = r_sections(k);
        r2 = r_sections(k+1);
        
        x1 = r1 * cos(angle_base);
        y1 = r1 * sin(angle_base);
        x2 = r2 * cos(angle_base);
        y2 = r2 * sin(angle_base);
        
        width = chord(k);
        dx = width * sin(angle_base);
        dy = -width * cos(angle_base);
        
        X = [x1, x2, x2+dx, x1+dx] * 1000;
        Y = [y1, y2, y2+dy, y1+dy] * 1000;
        
        stress_val = stress_centrifugal(j,k);
        stress_normalized = (stress_val - min(stress_centrifugal(:))) / (max(stress_centrifugal(:)) - min(stress_centrifugal(:)));
        color_idx = max(1, min(256, round(stress_normalized * 255 + 1)));
        face_color = stress_colormap(color_idx, :);
        
        patch(X, Y, face_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
end

draw_circle(0, 0, R_hub*1000, 'k', 'LineWidth', 2);
draw_circle(0, 0, R_max*1000, 'k--', 'LineWidth', 1);

xlabel('X (mm)');
ylabel('Y (mm)');
title('Centrifugal Stress Distribution');
colormap(jet);
c = colorbar;
c.Label.String = 'Centrifugal Stress (MPa)';
caxis([min(stress_centrifugal(:))/1e6, max(stress_centrifugal(:))/1e6]);
axis equal;
grid on;

% Figure 4: Bending Stress
figure('Name', 'Bending Stress Distribution', 'NumberTitle', 'off', 'Position', [200 200 1000 900]);

for j = 1:n_spokes
    angle_base = (j-1) * 2*pi/n_spokes;
    
    for k = 1:n_sections-1
        r1 = r_sections(k);
        r2 = r_sections(k+1);
        
        x1 = r1 * cos(angle_base);
        y1 = r1 * sin(angle_base);
        x2 = r2 * cos(angle_base);
        y2 = r2 * sin(angle_base);
        
        width = chord(k);
        dx = width * sin(angle_base);
        dy = -width * cos(angle_base);
        
        X = [x1, x2, x2+dx, x1+dx] * 1000;
        Y = [y1, y2, y2+dy, y1+dy] * 1000;
        
        stress_val = stress_bending(j,k);
        if max(stress_bending(:)) > 0
            stress_normalized = (stress_val - min(stress_bending(:))) / (max(stress_bending(:)) - min(stress_bending(:)));
        else
            stress_normalized = 0;
        end
        color_idx = max(1, min(256, round(stress_normalized * 255 + 1)));
        face_color = stress_colormap(color_idx, :);
        
        patch(X, Y, face_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
end

draw_circle(0, 0, R_hub*1000, 'k', 'LineWidth', 2);
draw_circle(0, 0, R_max*1000, 'k--', 'LineWidth', 1);

xlabel('X (mm)');
ylabel('Y (mm)');
title('Bending Stress Distribution');
colormap(jet);
c = colorbar;
c.Label.String = 'Bending Stress (MPa)';
caxis([min(stress_bending(:))/1e6, max(stress_bending(:))/1e6]);
axis equal;
grid on;

% Figure 5: Safety Factor
figure('Name', 'Safety Factor Distribution', 'NumberTitle', 'off', 'Position', [250 250 1000 900]);

for j = 1:n_spokes
    angle_base = (j-1) * 2*pi/n_spokes;
    
    for k = 1:n_sections-1
        r1 = r_sections(k);
        r2 = r_sections(k+1);
        
        x1 = r1 * cos(angle_base);
        y1 = r1 * sin(angle_base);
        x2 = r2 * cos(angle_base);
        y2 = r2 * sin(angle_base);
        
        width = chord(k);
        dx = width * sin(angle_base);
        dy = -width * cos(angle_base);
        
        X = [x1, x2, x2+dx, x1+dx] * 1000;
        Y = [y1, y2, y2+dy, y1+dy] * 1000;
        
        SF_val = safety_factor(j,k);
        
        if SF_val < 1.0
            face_color = [1 0 0];
        elseif SF_val < safety_factor_target
            face_color = [1 1 0];
        elseif SF_val < 4.0
            face_color = [0 1 0];
        else
            face_color = [0 0.5 1];
        end
        
        patch(X, Y, face_color, 'EdgeColor', 'k', 'LineWidth', 0.5);
        hold on;
    end
end

draw_circle(0, 0, R_hub*1000, 'k', 'LineWidth', 2);
draw_circle(0, 0, R_max*1000, 'k--', 'LineWidth', 1);

xlabel('X (mm)');
ylabel('Y (mm)');
title('Safety Factor Distribution');
axis equal;
grid on;

legend_x = R_max * 1.3 * 1000;
legend_y_start = R_max * 0.8 * 1000;
legend_spacing = 15;

rectangle('Position', [legend_x-5, legend_y_start-5, 10, 10], 'FaceColor', [1 0 0], 'EdgeColor', 'k');
text(legend_x+10, legend_y_start, 'SF < 1.0 (Critical)', 'FontSize', 10);

rectangle('Position', [legend_x-5, legend_y_start-legend_spacing-5, 10, 10], 'FaceColor', [1 1 0], 'EdgeColor', 'k');
text(legend_x+10, legend_y_start-legend_spacing, sprintf('SF < %.1f (Warning)', safety_factor_target), 'FontSize', 10);

rectangle('Position', [legend_x-5, legend_y_start-2*legend_spacing-5, 10, 10], 'FaceColor', [0 1 0], 'EdgeColor', 'k');
text(legend_x+10, legend_y_start-2*legend_spacing, 'SF < 4.0 (Good)', 'FontSize', 10);

rectangle('Position', [legend_x-5, legend_y_start-3*legend_spacing-5, 10, 10], 'FaceColor', [0 0.5 1], 'EdgeColor', 'k');
text(legend_x+10, legend_y_start-3*legend_spacing, 'SF > 4.0 (Excellent)', 'FontSize', 10);

% Figure 6: Single Blade Detail
figure('Name', 'Single Blade Stress Profile', 'NumberTitle', 'off', 'Position', [300 300 1200 600]);

subplot(1,2,1);
hold on;
for k = 1:n_sections-1
    r1 = r_sections(k);
    r2 = r_sections(k+1);
    
    stress_val = stress_total(1,k);
    stress_normalized = (stress_val - min(stress_total(:))) / (max(stress_total(:)) - min(stress_total(:)));
    color_idx = max(1, min(256, round(stress_normalized * 255 + 1)));
    face_color = stress_colormap(color_idx, :);
    
    rectangle('Position', [r1*1000, 0, (r2-r1)*1000, chord(k)*1000], ...
              'FaceColor', face_color, 'EdgeColor', 'k', 'LineWidth', 1);
end
xlabel('Radial Position (mm)');
ylabel('Chord (mm)');
title('Single Blade - Stress Distribution');
colormap(jet);
colorbar;
caxis([min(stress_total(:))/1e6, max(stress_total(:))/1e6]);
grid on;

subplot(1,2,2);
plot(r_sections*1000, stress_total(1,:)/1e6, 'r-', 'LineWidth', 3);
hold on;
plot(r_sections*1000, stress_centrifugal(1,:)/1e6, 'b--', 'LineWidth', 2);
plot(r_sections*1000, stress_bending(1,:)/1e6, 'g--', 'LineWidth', 2);
yline(sigma_yield/1e6, 'k--', 'Yield', 'LineWidth', 2);
xlabel('Radial Position (mm)');
ylabel('Stress (MPa)');
title('Stress Components Along Blade');
legend('Total', 'Centrifugal', 'Bending', 'Yield', 'Location', 'best');
grid on;

% Figure 7: Stress Contour Map
figure('Name', 'Stress Contour Map', 'NumberTitle', 'off', 'Position', [50 50 1400 600]);

stress_matrix = stress_total / 1e6;
[R_grid, Spoke_grid] = meshgrid(r_sections*1000, 1:n_spokes);

contourf(R_grid, Spoke_grid, stress_matrix, 20, 'LineColor', 'none');
colormap(jet);
c = colorbar;
c.Label.String = 'Total Stress (MPa)';
xlabel('Radial Position (mm)');
ylabel('Blade Number');
title('Stress Contour Map - All Blades');
hold on;
contour(R_grid, Spoke_grid, stress_matrix, [sigma_yield/1e6, sigma_yield/1e6], ...
        'LineColor', 'r', 'LineWidth', 3, 'LineStyle', '--');
grid on;

fprintf('\n✅ Stress visualization complete!\n');

fprintf('Generated 7 detailed stress distribution figures\n');
