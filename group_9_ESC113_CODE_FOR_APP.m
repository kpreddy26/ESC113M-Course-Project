classdef group_9_ESC113_CODE_FOR_APP < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        BUBBLEPOINTDEWPOINTCURVEButton  matlab.ui.control.Button
        CLEARPLOTButton                 matlab.ui.control.Button
        InitialTemperatureGuess90CSpinner  matlab.ui.control.Spinner
        InitialTemperatureGuess90CSpinnerLabel  matlab.ui.control.Label
        GROUP9Label                     matlab.ui.control.Label
        ToleranceValueEditField         matlab.ui.control.NumericEditField
        ToleranceValueEditFieldLabel    matlab.ui.control.Label
        PartialPressure2mmofHgEditField  matlab.ui.control.NumericEditField
        PartialPressure2mmofHgEditFieldLabel  matlab.ui.control.Label
        PartialPressure1mmofHgEditField  matlab.ui.control.NumericEditField
        PartialPressure1mmHgLabel       matlab.ui.control.Label
        Chemical2DropDown               matlab.ui.control.DropDown
        Chemical2DropDownLabel          matlab.ui.control.Label
        Chemical1DropDown               matlab.ui.control.DropDown
        Chemical1DropDownLabel          matlab.ui.control.Label
        RESULTSTextArea                 matlab.ui.control.TextArea
        RESULTSTextAreaLabel            matlab.ui.control.Label
        ApproximationTypeButtonGroup    matlab.ui.container.ButtonGroup
        NewtonRaphsonButton             matlab.ui.control.RadioButton
        RegulaFalsiButton               matlab.ui.control.RadioButton
        BisectionMethodButton           matlab.ui.control.RadioButton
        NewtonsMethodButton             matlab.ui.control.RadioButton
        PLOTButton                      matlab.ui.control.Button
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        coeffs1=zeros(1,3);
        coeffs2=zeros(1,3);
        Chemical1HasChanged = false; 
        Chemical2HasChanged = false;
    end
   
    methods (Access = private)
        
        function coeffs = getAntoineConstants(app, chemical)
            %values for the pre-fed chemicals
            switch chemical
                case 'Ethylbenzene'
                    coeffs = [6.95719, 1424.255, 213.21];
                case 'Toluene'
                    coeffs = [6.95464, 1344.8, 219.48];
                case 'Benzene'
                    coeffs = [6.90565, 1211.033, 220.79];
                case 'Acetone'
                    coeffs = [7.02447, 1161, 224];
                otherwise
                    coeffs = [NaN, NaN, NaN];
            end
        end
        
       
       
        %% Newton-Raphson
        function [x1, x2, T, errors] = estimateCompositionAndTemperature(app, A1, B1, C1, A2, B2, C2, P1, P2)
            % Newton-Raphson to solve for x1 and T 
            T0 = app.InitialTemperatureGuess90CSpinner.Value;
            x_guess = P1 / (P1 + P2);
            tol = app.ToleranceValueEditField.Value;  %tolerance input
            imax = 100;  %maximum number of iterations
            P = P1 + P2;  %total pressure
            y1 = P1 / P;  %vapour mole fractions
            y2 = P2 / P;
            %combined vector
            X0 = [x_guess; T0];
            errors = []; %to store erros for plotting

            for i = 1:imax
                %Antoines equations
                P_sat1 = 10^(A1 - B1 / (X0(2) + C1));
                P_sat2 = 10^(A2 - B2 / (X0(2) + C2));
                %Residual vector
                R = [X0(1) * P_sat1 - y1 * P;
                    (1 - X0(1)) * P_sat2 - y2 * P];
                %Derivatives for Jacobian
                dP1_dT = log(10) * P_sat1 * B1 / (X0(2) + C1)^2;
                dP2_dT = log(10) * P_sat2 * B2 / (X0(2) + C2)^2;
                %Jacobian analytically calculated
                J = [P_sat1, X0(1) * dP1_dT;
                    -P_sat2, (1 - X0(1)) * dP2_dT];
                %delta calculation
                del = -J \ R;
                X = X0 + del;  %updating X
                err = norm(X - X0, 2);
                errors(end+1) = err;
                %convergence criterion
                if err < tol
                    break;
                end
                X0 = X;
           end
            %storing final values
           x1 = X0(1);
           x2 = 1 - x1;
           T = X0(2);
        end
        
        %% Newton's Method
        function [T_newton, errors_newton] = NewtonMethod(app, f_T, df_T, imax, tol)
            T_newton = app.InitialTemperatureGuess90CSpinner.Value;  %initial guess
            errors_newton = [];  %storing error
            
            for i = 1:imax
                % calculation of f and its derivative for newton's method
                f_val = f_T(T_newton);
                df_val = df_T(T_newton);
                %prevent from dividing from near zero derivative
                if abs(df_val) < 1E-10
                    break;
                end
                %updating T
                T_newton = T_newton - f_val / df_val;
                err = abs(f_val);
                errors_newton = [errors_newton, err]; %storing errors for plotting
                %convergence criterion
                if err < tol
                    break;
                end
            end
        end
        
        %% **Bisection Method**
        function [T_bisection, errors_bisection] = BisectionMethod(app, f_T, tol, imax)
            %initializing values of bracket using both sides of initial
            %temperature guess
            T_low = app.InitialTemperatureGuess90CSpinner.Value - 40;
            T_high = app.InitialTemperatureGuess90CSpinner.Value + 60;
            errors_bisection = [];
            
            for i = 1:imax
                %convergence criterion
                if abs((T_high - T_low) / 2) <= tol
                    break;
                end
                %temperature calculation
                T_mid = (T_low + T_high) / 2;
                f_mid = f_T(T_mid);
                
                err = abs(f_mid);
                errors_bisection = [errors_bisection, err]; %storing errors for plotting
                %checking conditions for bisection method whether the root
                %lies in that bracket o/w shifting the bracket
                if f_T(T_low) * f_mid < 0
                    T_high = T_mid;
                else
                    T_low = T_mid;
                end
            end
            T_bisection = (T_low + T_high) / 2;
        end
        
        %% **Regula Falsi Method**
        function [T_regula, errors_regula] = RegulaFalsiMethod(app, f_T, tol, imax)
            %initializing values of bracket using both sides of initial
            %temperature guess
            T_low = app.InitialTemperatureGuess90CSpinner.Value - 40;
            T_high = app.InitialTemperatureGuess90CSpinner.Value + 60;
            errors_regula = [];
            
            for i = 1:imax
                f_low = f_T(T_low);
                f_high = f_T(T_high);
                %calculation of temperature using regula falsi formula
                T_regula = T_high - f_high * (T_high - T_low) / (f_high - f_low);
                f_regula = f_T(T_regula);
                
                err = abs(f_regula);
                errors_regula = [errors_regula, err]; %storing errors for plotting
                %convergence checking and conditions for bracket updation
                %of regula falsi
                if abs(f_regula) < tol
                    break;
                elseif f_T(T_low) * f_regula < 0
                    T_high = T_regula;
                else
                    T_low = T_regula;
                end
            end
        end
    
    %% Bubble Point calculations
    function T = solveBubblePointTemperature(app, A1, B1, C1, A2, B2, C2, x1, P)
    % Solve for T using Newton-Raphson: x1*Psat1 + x2*Psat2 = P
    tol = app.ToleranceValueEditField.Value;
    imax = 100;
    x2 = 1 - x1;  %sum of x1+x2=1
    T0 = app.InitialTemperatureGuess90CSpinner.Value;  %using initial temperature guess
    T = T0;

    for i = 1:imax
        %using antoines equations
        Psat1 = 10^(A1 - B1 / (T + C1));
        Psat2 = 10^(A2 - B2 / (T + C2));
        %residual equation
        f = x1 * Psat1 + x2 * Psat2 - P;
        %derivatives for newton's method
        dPsat1_dT = log(10) * Psat1 * B1 / (T + C1)^2;
        dPsat2_dT = log(10) * Psat2 * B2 / (T + C2)^2;
        df = x1 * dPsat1_dT + x2 * dPsat2_dT;
        %updating T
        delta = -f / df;
        T = T + delta;
        %convergence criterion
        if abs(delta) < tol
            break;
        end
    end
    end

    %% Dew pointt calculations
    function T = solveDewPointTemperature(app, A1, B1, C1, A2, B2, C2, y1, P)
    tol = app.ToleranceValueEditField.Value;
    imax = 100;
    y2 = 1 - y1;  %mole fraction sum=1
    T = app.InitialTemperatureGuess90CSpinner.Value;

    for i = 1:imax
        %using antoines equation
        Psat1 = 10^(A1 - B1 / (T + C1));
        Psat2 = 10^(A2 - B2 / (T + C2));
        %residual function
        f = y1 * P / Psat1 + y2 * P / Psat2 - 1;
        %derivatives for newton's method
        dPsat1_dT = log(10) * Psat1 * B1 / (T + C1)^2;
        dPsat2_dT = log(10) * Psat2 * B2 / (T + C2)^2;
        df = -y1 * P * dPsat1_dT / Psat1^2 - y2 * P * dPsat2_dT / Psat2^2;
        %updating T
        delta = -f / df;
        T = T + delta;
        %convergence criterion
        if abs(delta) < tol
            break;
        end
    end
    end


        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: PLOTButton
        function PLOTButtonPushed(app, event)
            if strcmp(app.Chemical1DropDown.Value, app.Chemical2DropDown.Value)
                uialert(app.UIFigure, 'Please select two different chemicals.', 'Input Error');
                return;
            end
            
            %% Given Partial Pressures
            P1 = app.PartialPressure1mmofHgEditField.Value;     % Partial Pressure of Chemical 1
            P2 = app.PartialPressure2mmofHgEditField.Value;     % Partial Pressure of Chemical 2
            P = (P1 + P2); % Total pressure assuming equal mole fractions
           
            % Antoine Coefficients
            % default conditions
            if ~app.Chemical1HasChanged
                chemical = app.Chemical1DropDown.Value;
                app.coeffs1 = app.getAntoineConstants(chemical);
            end

            if ~app.Chemical2HasChanged
                chemical = app.Chemical2DropDown.Value;
                app.coeffs2 = app.getAntoineConstants(chemical);
            end
            %Assigning values to coefficients
            A1 = app.coeffs1(1);
            B1 = app.coeffs1(2);
            C1 = app.coeffs1(3);

            A2 = app.coeffs2(1);
            B2 = app.coeffs2(2);
            C2 = app.coeffs2(3);
            %tolerance value input from user
            tol = app.ToleranceValueEditField.Value;
            imax = 100;  %maximum number of iterations
            T=app.InitialTemperatureGuess90CSpinner.Value;  %initial T guess from user
            %% Estimation of x1 using Newton-Raphson**
            [x1, x2, T0, errors] = app.estimateCompositionAndTemperature(A1, B1, C1, A2, B2, C2, P1, P2);

            %% Estimation of Temperature Using Newton, Bisection, and Regula Falsi**
            
            %Function of single variable T for further T estimation using
            %shorthand notation for function definition
            f_T = @(T) x1 * 10^(A1 - B1 / (T + C1)) + x2 * 10^(A2 - B2 / (T + C2)) - P;
            df_T = @(T) log(10) * (x1 * 10^(A1 - B1 / (T + C1)) * B1 / (T + C1)^2 + ...
                        x2 * 10^(A2 - B2 / (T + C2)) * B2 / (T + C2)^2);

            %% **Newton's Method**
            [T_newton, errors_newton] = NewtonMethod(app, f_T, df_T, imax, tol);

            %% **Bisection Method**
            [T_bisection, errors_bisection] = BisectionMethod(app, f_T, tol, imax);

            %% **Regula Falsi**
            [T_regula, errors_regula] = RegulaFalsiMethod(app, f_T, tol, imax);

            %% **Convergence Plots
            %plotting according to radio button selection
            if app.NewtonsMethodButton.Value == true
                x=1:length(errors_newton);
                y=errors_newton;
                plot(app.UIAxes, x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName=sprintf("Newton's Method @ T=%.1f°C", T));
            elseif app.BisectionMethodButton.Value == true
                x=1:length(errors_bisection);
                y=errors_bisection;
                plot(app.UIAxes, x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName=sprintf("Bisection Method @ T=%.1f°C", T));
            elseif app.RegulaFalsiButton.Value == true
                x=1:length(errors_regula);
                y=errors_regula;
                plot(app.UIAxes, x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName=sprintf("Regula Falsi @ T=%.1f°C", T));
            elseif app.NewtonRaphsonButton.Value == true
                x=1:length(errors);
                y=errors;
                plot(app.UIAxes, x, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, DisplayName=sprintf("Newton Raphson @ T=%.1f°C", T));
            end
            %features of the plot
            legend(app.UIAxes)
            grid(app.UIAxes, 'on');
            hold(app.UIAxes,'on');
            title(app.UIAxes, 'Error Convergence Plot');
            xlabel(app.UIAxes, 'Iteration');
            ylabel(app.UIAxes, 'Absolute Error');
            %result string to be displayed in the RESULTS area
            resultStr = sprintf(['Using Newton-Raphson:\nEstimated x1 (%s): %.4f\nEstimated x2 (%s): %.4f\nEstimated T: %.4f °C\n\n' ...
                     'Newtons Method: Temperature = %.4f °C\nBisection Method: Temperature = %.4f °C\nRegula Falsi Method: Temperature = %.4f °C'], ...
                     app.Chemical1DropDown.Value,x1,app.Chemical2DropDown.Value, x2, T0, T_newton, T_bisection, T_regula);
            app.RESULTSTextArea.Value = splitlines(resultStr);
        end

        % Value changed function: Chemical1DropDown
        function Chemical1DropDownValueChanged(app, event)
            %to update the value of chemical 1's Antoines coefficients if changed
            chemical = app.Chemical1DropDown.Value;
            coeffs = app.getAntoineConstants(chemical);
            app.coeffs1 = coeffs;
            app.Chemical1HasChanged = true; %flag if changed
        end

        % Value changed function: Chemical2DropDown
        function Chemical2DropDownValueChanged(app, event)
            %to update the value of chemical 2's Antoines coefficients if changed
            chemical = app.Chemical2DropDown.Value;
            coeffs = app.getAntoineConstants(chemical);
            app.coeffs2 = coeffs;
            app.Chemical2HasChanged = true;  %flag if changed
        end

        % Value changed function: InitialTemperatureGuess90CSpinner
        function InitialTemperatureGuess90CSpinnerValueChanged(app, event)
            %to update the value of T as per user selection
            value = app.InitialTemperatureGuess90CSpinner.Value;
            app.InitialTemperatureGuess90CSpinnerLabel.Text = ...
            sprintf('Initial Temperature Guess: \n%.1f °C', value);
        end

        % Button pushed function: CLEARPLOTButton
        function CLEARPLOTButtonPushed(app, event)
            %to clear the axes of plot
            cla(app.UIAxes);
        end

        % Button pushed function: BUBBLEPOINTDEWPOINTCURVEButton
        function BUBBLEPOINTDEWPOINTCURVEButtonPushed(app, event)
            if strcmp(app.Chemical1DropDown.Value, app.Chemical2DropDown.Value)
                uialert(app.UIFigure, 'Please select two different chemicals.', 'Input Error');
                return;
            end
            % Given Partial Pressures
            P1 = app.PartialPressure1mmofHgEditField.Value;     % Partial Pressure of Chemical 1
            P2 = app.PartialPressure2mmofHgEditField.Value;     % Partial Pressure of Chemical 2
            P = (P1 + P2); % Total pressure assuming equal mole fractions
            if ~app.Chemical1HasChanged
                chemical = app.Chemical1DropDown.Value;
                app.coeffs1 = app.getAntoineConstants(chemical);
            end

            if ~app.Chemical2HasChanged
                chemical = app.Chemical2DropDown.Value;
                app.coeffs2 = app.getAntoineConstants(chemical);
            end
            %Assigning values to coefficients
            A1 = app.coeffs1(1);
            B1 = app.coeffs1(2);
            C1 = app.coeffs1(3);

            A2 = app.coeffs2(1);
            B2 = app.coeffs2(2);
            C2 = app.coeffs2(3);
            x1 = linspace(0, 1, 50);
            T_bubble = zeros(size(x1));
            %calculating T bubble for x1 ranging from 0 to 1
            for i = 1:length(x1)
                x0 = x1(i);
                T_bubble(i)=app.solveBubblePointTemperature(A1, B1, C1, A2, B2, C2, x0, P);
            end

            y1 = linspace(0, 1, 50);
            T_dew = zeros(size(y1));
  
            %calculating T dew for y1 changing from 0 to 1
            for i = 1:length(y1)
                y0 = y1(i);
                T_dew(i) = app.solveDewPointTemperature(A1, B1, C1, A2, B2, C2, y0, P);
            end
            %getting the calculated temperature of the mixture displayed
            %earlier
            [~, ~, T_isotherm, ~] = app.estimateCompositionAndTemperature(A1, B1, C1, A2, B2, C2, P1, P2);
            % Plotting the dew point and bubble point curves
            plot(app.UIAxes2, x1, T_bubble, 'b-', 'DisplayName', 'Bubble Point');
            hold(app.UIAxes2, 'on');
            plot(app.UIAxes2, y1, T_dew, 'r-', 'DisplayName', 'Dew Point');
            % Plotting the Isotherm
            plot(app.UIAxes2,[0 1],[T_isotherm T_isotherm], 'g-', 'DisplayName', 'Isotherm');
            %plot features
            legend(app.UIAxes2, 'show');
            xlabel(app.UIAxes2, 'Mole Fraction');
            ylabel(app.UIAxes2, 'Temperature (°C)');
            title(app.UIAxes2, 'Bubble & Dew Point Curves');
            grid(app.UIAxes2, 'on');
            hold(app.UIAxes2, 'off');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.7373 0.8353 0.8392];
            app.UIFigure.Position = [100 100 1013 671];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Error Convergence Plot')
            xlabel(app.UIAxes, 'Iterations')
            ylabel(app.UIAxes, 'Error')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontName = 'Times New Roman';
            app.UIAxes.FontWeight = 'bold';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.XMinorGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.YMinorGrid = 'on';
            app.UIAxes.FontSize = 14;
            app.UIAxes.Position = [629 348 338 250];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'Bubble Point - Dew Point Curve')
            xlabel(app.UIAxes2, 'Mole Fraction')
            ylabel(app.UIAxes2, 'Temperature (°C)')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.FontName = 'Times New Roman';
            app.UIAxes2.FontWeight = 'bold';
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.XMinorGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.YMinorGrid = 'on';
            app.UIAxes2.Position = [17 60 343 248];

            % Create PLOTButton
            app.PLOTButton = uibutton(app.UIFigure, 'push');
            app.PLOTButton.ButtonPushedFcn = createCallbackFcn(app, @PLOTButtonPushed, true);
            app.PLOTButton.FontName = 'Times New Roman';
            app.PLOTButton.FontSize = 14;
            app.PLOTButton.FontWeight = 'bold';
            app.PLOTButton.Position = [348 396 98 26];
            app.PLOTButton.Text = 'PLOT';

            % Create ApproximationTypeButtonGroup
            app.ApproximationTypeButtonGroup = uibuttongroup(app.UIFigure);
            app.ApproximationTypeButtonGroup.TitlePosition = 'centertop';
            app.ApproximationTypeButtonGroup.Title = 'Approximation Type';
            app.ApproximationTypeButtonGroup.BackgroundColor = [0.949 0.8157 0.8157];
            app.ApproximationTypeButtonGroup.FontName = 'Times New Roman';
            app.ApproximationTypeButtonGroup.FontWeight = 'bold';
            app.ApproximationTypeButtonGroup.FontSize = 14;
            app.ApproximationTypeButtonGroup.Position = [345 447 249 114];

            % Create NewtonsMethodButton
            app.NewtonsMethodButton = uiradiobutton(app.ApproximationTypeButtonGroup);
            app.NewtonsMethodButton.Text = 'Newton''s Method';
            app.NewtonsMethodButton.FontName = 'Times New Roman';
            app.NewtonsMethodButton.Position = [11 66 107 22];
            app.NewtonsMethodButton.Value = true;

            % Create BisectionMethodButton
            app.BisectionMethodButton = uiradiobutton(app.ApproximationTypeButtonGroup);
            app.BisectionMethodButton.Text = 'Bisection Method';
            app.BisectionMethodButton.FontName = 'Times New Roman';
            app.BisectionMethodButton.Position = [11 44 107 22];

            % Create RegulaFalsiButton
            app.RegulaFalsiButton = uiradiobutton(app.ApproximationTypeButtonGroup);
            app.RegulaFalsiButton.Text = 'Regula Falsi';
            app.RegulaFalsiButton.FontName = 'Times New Roman';
            app.RegulaFalsiButton.Position = [11 22 82 22];

            % Create NewtonRaphsonButton
            app.NewtonRaphsonButton = uiradiobutton(app.ApproximationTypeButtonGroup);
            app.NewtonRaphsonButton.Text = 'Newton Raphson';
            app.NewtonRaphsonButton.FontName = 'Times New Roman';
            app.NewtonRaphsonButton.Position = [11 0 105 22];

            % Create RESULTSTextAreaLabel
            app.RESULTSTextAreaLabel = uilabel(app.UIFigure);
            app.RESULTSTextAreaLabel.HorizontalAlignment = 'center';
            app.RESULTSTextAreaLabel.FontName = 'Times New Roman';
            app.RESULTSTextAreaLabel.FontSize = 14;
            app.RESULTSTextAreaLabel.FontWeight = 'bold';
            app.RESULTSTextAreaLabel.Position = [570 255 198 30];
            app.RESULTSTextAreaLabel.Text = 'RESULTS';

            % Create RESULTSTextArea
            app.RESULTSTextArea = uitextarea(app.UIFigure);
            app.RESULTSTextArea.FontName = 'Times New Roman';
            app.RESULTSTextArea.FontSize = 14;
            app.RESULTSTextArea.FontWeight = 'bold';
            app.RESULTSTextArea.FontAngle = 'italic';
            app.RESULTSTextArea.Position = [392 94 554 162];

            % Create Chemical1DropDownLabel
            app.Chemical1DropDownLabel = uilabel(app.UIFigure);
            app.Chemical1DropDownLabel.FontName = 'Times New Roman';
            app.Chemical1DropDownLabel.FontSize = 14;
            app.Chemical1DropDownLabel.FontWeight = 'bold';
            app.Chemical1DropDownLabel.Position = [28 539 72 22];
            app.Chemical1DropDownLabel.Text = 'Chemical 1';

            % Create Chemical1DropDown
            app.Chemical1DropDown = uidropdown(app.UIFigure);
            app.Chemical1DropDown.Items = {'Ethylbenzene', 'Toluene', 'Benzene', 'Acetone'};
            app.Chemical1DropDown.ValueChangedFcn = createCallbackFcn(app, @Chemical1DropDownValueChanged, true);
            app.Chemical1DropDown.Position = [133 539 167 22];
            app.Chemical1DropDown.Value = 'Ethylbenzene';

            % Create Chemical2DropDownLabel
            app.Chemical2DropDownLabel = uilabel(app.UIFigure);
            app.Chemical2DropDownLabel.FontName = 'Times New Roman';
            app.Chemical2DropDownLabel.FontSize = 14;
            app.Chemical2DropDownLabel.FontWeight = 'bold';
            app.Chemical2DropDownLabel.Position = [30 469 72 22];
            app.Chemical2DropDownLabel.Text = 'Chemical 2';

            % Create Chemical2DropDown
            app.Chemical2DropDown = uidropdown(app.UIFigure);
            app.Chemical2DropDown.Items = {'Toluene', 'Ethylbenzene', 'Acetone', 'Benzene'};
            app.Chemical2DropDown.ValueChangedFcn = createCallbackFcn(app, @Chemical2DropDownValueChanged, true);
            app.Chemical2DropDown.Position = [135 469 165 22];
            app.Chemical2DropDown.Value = 'Toluene';

            % Create PartialPressure1mmHgLabel
            app.PartialPressure1mmHgLabel = uilabel(app.UIFigure);
            app.PartialPressure1mmHgLabel.HorizontalAlignment = 'center';
            app.PartialPressure1mmHgLabel.VerticalAlignment = 'top';
            app.PartialPressure1mmHgLabel.FontName = 'Times New Roman';
            app.PartialPressure1mmHgLabel.FontSize = 14;
            app.PartialPressure1mmHgLabel.FontWeight = 'bold';
            app.PartialPressure1mmHgLabel.Position = [28 495 112 34];
            app.PartialPressure1mmHgLabel.Text = {'Partial Pressure 1'; '(mm of Hg)'};

            % Create PartialPressure1mmofHgEditField
            app.PartialPressure1mmofHgEditField = uieditfield(app.UIFigure, 'numeric');
            app.PartialPressure1mmofHgEditField.Limits = [0 1000];
            app.PartialPressure1mmofHgEditField.HorizontalAlignment = 'left';
            app.PartialPressure1mmofHgEditField.Position = [184 507 116 22];
            app.PartialPressure1mmofHgEditField.Value = 250;

            % Create PartialPressure2mmofHgEditFieldLabel
            app.PartialPressure2mmofHgEditFieldLabel = uilabel(app.UIFigure);
            app.PartialPressure2mmofHgEditFieldLabel.HorizontalAlignment = 'center';
            app.PartialPressure2mmofHgEditFieldLabel.FontName = 'Times New Roman';
            app.PartialPressure2mmofHgEditFieldLabel.FontSize = 14;
            app.PartialPressure2mmofHgEditFieldLabel.FontWeight = 'bold';
            app.PartialPressure2mmofHgEditFieldLabel.Position = [29 421 112 34];
            app.PartialPressure2mmofHgEditFieldLabel.Text = {'Partial Pressure 2'; '(mm of Hg)'};

            % Create PartialPressure2mmofHgEditField
            app.PartialPressure2mmofHgEditField = uieditfield(app.UIFigure, 'numeric');
            app.PartialPressure2mmofHgEditField.Limits = [0 1000];
            app.PartialPressure2mmofHgEditField.HorizontalAlignment = 'left';
            app.PartialPressure2mmofHgEditField.Position = [184 433 116 22];
            app.PartialPressure2mmofHgEditField.Value = 343;

            % Create ToleranceValueEditFieldLabel
            app.ToleranceValueEditFieldLabel = uilabel(app.UIFigure);
            app.ToleranceValueEditFieldLabel.FontName = 'Times New Roman';
            app.ToleranceValueEditFieldLabel.FontSize = 14;
            app.ToleranceValueEditFieldLabel.FontWeight = 'bold';
            app.ToleranceValueEditFieldLabel.Position = [28 391 101 22];
            app.ToleranceValueEditFieldLabel.Text = 'Tolerance Value';

            % Create ToleranceValueEditField
            app.ToleranceValueEditField = uieditfield(app.UIFigure, 'numeric');
            app.ToleranceValueEditField.Limits = [1e-12 0.1];
            app.ToleranceValueEditField.HorizontalAlignment = 'left';
            app.ToleranceValueEditField.Position = [184 391 116 22];
            app.ToleranceValueEditField.Value = 1e-06;

            % Create GROUP9Label
            app.GROUP9Label = uilabel(app.UIFigure);
            app.GROUP9Label.HorizontalAlignment = 'center';
            app.GROUP9Label.FontName = 'Times New Roman';
            app.GROUP9Label.FontSize = 20;
            app.GROUP9Label.FontWeight = 'bold';
            app.GROUP9Label.Position = [-22 573 989 112];
            app.GROUP9Label.Text = 'Computation of Liquid-Phase Composition and Temperature of a Vapour-Liquid Mixture';

            % Create InitialTemperatureGuess90CSpinnerLabel
            app.InitialTemperatureGuess90CSpinnerLabel = uilabel(app.UIFigure);
            app.InitialTemperatureGuess90CSpinnerLabel.HorizontalAlignment = 'center';
            app.InitialTemperatureGuess90CSpinnerLabel.FontName = 'Times New Roman';
            app.InitialTemperatureGuess90CSpinnerLabel.FontSize = 14;
            app.InitialTemperatureGuess90CSpinnerLabel.FontWeight = 'bold';
            app.InitialTemperatureGuess90CSpinnerLabel.Position = [345 316 196 34];
            app.InitialTemperatureGuess90CSpinnerLabel.Text = {'Initial Temperature Guess '; '90 °C'};

            % Create InitialTemperatureGuess90CSpinner
            app.InitialTemperatureGuess90CSpinner = uispinner(app.UIFigure);
            app.InitialTemperatureGuess90CSpinner.Limits = [30 170];
            app.InitialTemperatureGuess90CSpinner.ValueChangedFcn = createCallbackFcn(app, @InitialTemperatureGuess90CSpinnerValueChanged, true);
            app.InitialTemperatureGuess90CSpinner.HorizontalAlignment = 'left';
            app.InitialTemperatureGuess90CSpinner.Position = [549 314 81 34];
            app.InitialTemperatureGuess90CSpinner.Value = 90;

            % Create CLEARPLOTButton
            app.CLEARPLOTButton = uibutton(app.UIFigure, 'push');
            app.CLEARPLOTButton.ButtonPushedFcn = createCallbackFcn(app, @CLEARPLOTButtonPushed, true);
            app.CLEARPLOTButton.FontName = 'Times New Roman';
            app.CLEARPLOTButton.FontSize = 14;
            app.CLEARPLOTButton.FontWeight = 'bold';
            app.CLEARPLOTButton.Position = [476 396 118 26];
            app.CLEARPLOTButton.Text = 'CLEAR PLOT';

            % Create BUBBLEPOINTDEWPOINTCURVEButton
            app.BUBBLEPOINTDEWPOINTCURVEButton = uibutton(app.UIFigure, 'push');
            app.BUBBLEPOINTDEWPOINTCURVEButton.ButtonPushedFcn = createCallbackFcn(app, @BUBBLEPOINTDEWPOINTCURVEButtonPushed, true);
            app.BUBBLEPOINTDEWPOINTCURVEButton.FontName = 'Times New Roman';
            app.BUBBLEPOINTDEWPOINTCURVEButton.FontSize = 14;
            app.BUBBLEPOINTDEWPOINTCURVEButton.FontWeight = 'bold';
            app.BUBBLEPOINTDEWPOINTCURVEButton.Position = [71 324 267 26];
            app.BUBBLEPOINTDEWPOINTCURVEButton.Text = 'BUBBLE POINT - DEW POINT CURVE';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = group_9_ESC113_CODE_FOR_APP

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end