import { logger } from "../util/logger";
import { PrismaClient, RefreshToken } from "@prisma/client";
import { SHA256 } from "../util/utils";
import { generateJWT, generateRT } from "../util/jwt.utils";

class UserModel {
    prisma: PrismaClient;

    constructor() {
        this.prisma = new PrismaClient();
    }

    async findByEmail(email: string) {
        return this.prisma.user.findUnique({
            where: { email: email },
        });
    }

    async findUserIdByEmail(email: string) {
        return (
            await this.prisma.user.findUnique({
                where: { email: email },
            })
        )?.id;
    }
    async create(email: string, password: string) {
        try {
            const hashed_pw = await SHA256(password);
            return await this.prisma.user.create({
                data: {
                    email: email,
                    password: hashed_pw,
                },
            });
        } catch (err) {
            throw {
                stack: "error in creating user, email: " + email,
                status: 500,
                message: "ошибка при создании пользователя",
            };
        }
    }
}

export const userModel = new UserModel();

class RefreshTokenModel {
    prisma: PrismaClient;
    constructor() {
        this.prisma = new PrismaClient();
    }

    async findByEmail(email: string) {
        return this.prisma.refreshToken.findUnique({
            where: { userid: await userModel.findUserIdByEmail(email) },
        });
    }

    async create(userid: bigint) {
        try {
            const rt = await this.prisma.refreshToken.create({
                data: {
                    userid: userid,
                    token: generateRT(),
                },
            });
            logger.info("refresh token for user created: " + userid);
            return rt;
        } catch (err) {
            throw {
                stack: "Error in creating RefreshToken, userid: " + userid,
                status: 500,
                message: "Ошибка при создании токена",
            };
        }
    }

    async updateByUserId(userid: bigint) {
        try {
            const rt = await this.prisma.refreshToken.findUnique({
                where: { userid: userid },
            });
            if (rt) {
                rt.token = generateRT();
                logger.info("refresh token updated for user: " + userid);
                return rt;
            } else {
                return this.create(userid);
            }
        } catch (err) {
            throw {
                stack: "error in updating RefreshToken, userid: " + userid,
                status: 500,
                message: "ошибка при обновлении токена",
            };
        }
    }
    
    async deleteByUserId(userid: bigint){
        try {
            const rt = await this.prisma.refreshToken.delete({
                where: { userid: userid },
            });
            if (rt) {
                rt.token = generateRT();
                logger.info("refresh token deleted for user: " + userid);
                return rt;
            } else {
                return this.create(userid);
            }
        } catch (err) {
            throw {
                stack: "error in deleting RefreshToken, userid: " + userid,
                status: 500,
                message: "ошибка при обновлении токена",
            };
        }
    }
}

export const refreshTokenModel = new RefreshTokenModel();
