function y = impulse_function(type, t, parameters)

    switch type
        case "gaussian"
            mu = parameters(1);
            sig = parameters(2);
            amplitude = parameters(3);
            y = amplitude * exp(-((((t)-mu).^2)/(2*sig.^2)));
    end

end
