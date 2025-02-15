import { Request, Response, NextFunction } from "express";
import { errorResponce } from "../util/utils";

export const errorHandler = (
    err: any,
    req: Request,
    res: Response,
    next: NextFunction,
): void => {
    errorResponce(err, res);
};
