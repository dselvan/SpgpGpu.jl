function minimize(X, f, len, P1, P2, P3, P4)
    RHO = 0.01
    SIG = 0.5
    INT = 0.1
    EXT = 3.0
    MAX = 20
    RATIO = 100

    # quick bodge for demo
    argstr = string(f, "(X,P1,P2,P3,P4)")
    if max(size(len) == 2)
        red = len[2]
        len = len[1]
    else
        red = 1
    end

    if len > 0
        S = "Linesearch"
    else
        S = "Function evaluation"
    end

    i = 0
    ls_failed = 0
    fX = []
    [f1, df1] = eval(Meta.parse(argstr))
    i = i + (len < 0)
    s = -df1
    d1 = transpose(-s) * s
    z1 = red / (1 - d1)

    while i < abs(len)
        i = i + (len > 0)

        X0 = X
        f0 = f1
        df0 = df1

        X = X + x1 * s
        [f2, df2] = eval(Meta.parse(argstr))
        i = i + (len < 0)
        d2 = transpose(df2) * s
        f3 = f1
        d3 = d1
        z3 = -z1
        if len > 0
            M = MAX
        else
            M = min(MAX, -len - i)
        end
        success = 0
        limit = -1
        while true
            while ((f2 > f1 + z1 * RHO * d1) || (d2 > -SIG * d1)) && (M > 0)
                limit = z1
                if f2 > f1
                    z2 = z3 - (0.5 * d3 * z3 * z3) / (d3 * z3 + f2 - f3)
                else
                    A = 6 * (f2 - f3) / z3 + 3 * (d2 + d3)
                    B = 3 * (f3 - f2) - z3 * (d3 + 2 * d2)
                    z2 = (sqrt(B * B - A * d2 * z3 * z3) - B) / A'
                end
                if isnan(z2) || isinf(z2)
                    z2 = z3 / 2
                end
                z2 = max(min(z2, INT * z3), (1 - INT) * z3)
                z1 = z1 + z2
                X = X + z2 * s
                [f2, df2] = eval(Meta.parse(argstr))
                M = M - 1
                i = i + (len < 0)
                d2 = transpose(df2) * s
                z3 = z3 - z2
            end
            if f2 > f1 + z1 * RHO * d1 || d2 > -SIG * d1
                break
            elseif d2 > SIG * d1
                success = true;
                break
            elseif M == 0
                break
            end
            A = 6 * (f2 - f3) / z3 + 3 * (d2 + d3)
            B = 3 * (f3 - f2) - z3 * (d3 + 2 * d2)
            z2 = -d2 * z3 * z3 / (B + sqrt(B * B - A * d2 * z3 * z3))
            if ~isreal(z2) | isnan(z2) | isinf(z2) | z2 < 0
                if limit < -0.5
                    z2 = z1 * (EXT - 1)
                    z2 = (limit - z1) / 2
                end
            elseif (limit > -0.5) && (z2 + z1 > limit)
                z2 = (limit - z1) / 2
            elseif (limit < -0.5) && (z2 + z1 > z1 * EXT)
                z2 = z1 * (EXT - 1.0)
            elseif z2 < -z3 * INT
                z2 = -z3 * INT
            elseif (limit > -0.5) && (z2 < (limit - z1) * (1.0 - INT))
                z2 = (limit - z1) * (1.0 - INT)
            end
            f3 = f2
            d3 = d2
            z3 = -z2
            z1 = z1 + z2
            X = X + z2 * s
            [f2, df2] = eval(Meta.parse(argstr))
            M = M - 1
            i = i + (length < 0)
            d2 = df2' * s
        end

        if success
            f1 = f2;
            fX = [transpose(fX) f1];
            @printf "%s %6i; Value %4.6e\n" S i f1
            s = (transpose(df2)*df2-transpose(df1)*df2)/(transpose(df1)*df1)*s - df2;
            tmp = df1;
            df1=df2;
            df2=tmp;
            d2 = transpose(df1)*s;
            if d2 > 0
                s = -df1;
                d2 = -transpose(s)*s;
            end
            z1 = z1 * min(RATIO, d1/(d2-typemin(Int64)));
            d1=d2;
            ls_failed = 0;
        else
            X = X0;
            f1 = f0;
            df1 = df0;
            if ls_failed || i > abs(len)
                break;
            end
            tmp = df1;
            df1 = df2;
            df2 = tmp;
            s = -df1;
            d1 = -transpose(s)*s;
            z1 = 1/(1-d1);
            ls_failed = 1;
        end       
    end
    @printf "\n"
    return X, fX, i
end
